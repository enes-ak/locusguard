"""Per-locus assignment orchestrator.

For each locus, fetch reads overlapping the primary + paralog regions,
build AnalyzedRead objects, cluster them by PSV-pattern, run evidence
adapters (including haplotype_consistency on clustered reads), score, and
emit Assignments tied back to their haplotype cluster.
"""
from __future__ import annotations

import hashlib
from collections.abc import Iterator

import pysam

from locusguard.cn import estimate_cn
from locusguard.config.resolver import ResolvedProfile, resolve_profile
from locusguard.config.schema import PSV, LocusConfig
from locusguard.depth import DepthStats, compute_region_depth
from locusguard.evidence.base import EvidenceSource, ReadTech
from locusguard.evidence.coverage_ratio import CoverageRatioEvidence
from locusguard.evidence.haplotype_consistency import HaplotypeConsistencyEvidence
from locusguard.evidence.mapq_pattern import MapqPatternEvidence
from locusguard.evidence.psv import PSVEvidence
from locusguard.evidence.softclip import SoftclipEvidence
from locusguard.haplotype import (
    assign_cluster_locus,
    cluster_reads,
    detect_gene_conversion,
)
from locusguard.io.bam import BamReader
from locusguard.io.fasta import FastaReader
from locusguard.scoring import score_assignment
from locusguard.types import (
    AnalyzedRead,
    Assignment,
    CnEstimate,
    EvidenceScore,
    HaplotypeCluster,
    PSVObs,
)

_CIGAR_SOFT_CLIP = 4


def build_analyzed_reads(
    bam: BamReader,
    fasta: FastaReader,
    config: LocusConfig,
    profile_name: str | None,
) -> Iterator[AnalyzedRead]:
    """Iterate reads overlapping primary + paralog regions.

    Only yield reads that physically reach at least one configured PSV
    position — reads without any PSV signal contribute nothing downstream
    and clutter outputs.
    """
    regions = [config.coordinates.primary, *config.coordinates.paralogs.values()]
    long_read_hint = bam.estimated_is_long_read()
    seen: set[str] = set()

    for region in regions:
        for read in bam.fetch(region.chrom, region.start, region.end):
            if read.query_name is None or read.query_name in seen:
                continue
            seen.add(read.query_name)
            analyzed = _analyze_read(read, config.psvs, long_read_hint)
            if any(obs.reach for obs in analyzed.psv_observations.values()):
                yield analyzed


def _analyze_read(
    read: pysam.AlignedSegment,
    psvs: list[PSV],
    long_read_hint: bool,
) -> AnalyzedRead:
    observations: dict[str, PSVObs] = {}
    ref_to_query = _ref_to_query_index(read)
    for psv in psvs:
        if psv.chrom != read.reference_name:
            observations[psv.name] = PSVObs(base="N", qual=0, reach=False)
            continue
        # Config uses 1-based genomics convention; pysam returns 0-based.
        q_idx = ref_to_query.get(psv.pos - 1)
        if q_idx is None:
            observations[psv.name] = PSVObs(base="N", qual=0, reach=False)
            continue
        base = read.query_sequence[q_idx].upper() if read.query_sequence else "N"
        qual = (
            read.query_qualities[q_idx]
            if read.query_qualities is not None and q_idx < len(read.query_qualities)
            else 0
        )
        observations[psv.name] = PSVObs(base=base, qual=int(qual), reach=True)

    softclip_5p, softclip_3p = _softclip_amounts(read)
    return AnalyzedRead(
        read_id=read.query_name or "",
        aligned_chrom=read.reference_name or "",
        aligned_pos=read.reference_start,
        psv_observations=observations,
        mapq=read.mapping_quality,
        softclip_5p=softclip_5p,
        softclip_3p=softclip_3p,
        is_long_read=long_read_hint,
        is_supplementary=read.is_supplementary,
        original_mapq_zero=read.mapping_quality == 0,
        cluster_consensus=None,
    )


def _ref_to_query_index(read: pysam.AlignedSegment) -> dict[int, int]:
    pairs = read.get_aligned_pairs(matches_only=True)
    return {ref: q for q, ref in pairs}


def _softclip_amounts(read: pysam.AlignedSegment) -> tuple[int, int]:
    cigar = read.cigartuples or []
    if not cigar:
        return 0, 0
    sc_5 = cigar[0][1] if cigar[0][0] == _CIGAR_SOFT_CLIP else 0
    sc_3 = cigar[-1][1] if cigar[-1][0] == _CIGAR_SOFT_CLIP else 0
    return sc_5, sc_3


class LocusAssigner:
    """Orchestrates clustering + evidence collection + scoring for one locus."""

    def __init__(self, config: LocusConfig, profile_name: str | None) -> None:
        self._config = config
        self._profile_name = profile_name
        self._profile: ResolvedProfile = resolve_profile(
            config.evidence_weights, profile_name
        )
        self._locus_key = self._compute_locus_key()
        self._psv_names = [p.name for p in config.psvs]
        self._haplotype_clusters: list[HaplotypeCluster] = []
        self._cn_estimate: CnEstimate | None = None
        # _adapters is rebuilt per-assign because CoverageRatioEvidence needs
        # a fresh depths_by_locus from the current BAM.
        self._base_adapters: list[EvidenceSource] = [
            PSVEvidence(target_locus=config.locus.id),
            MapqPatternEvidence(),
            SoftclipEvidence(),
            HaplotypeConsistencyEvidence(),
        ]

    @property
    def warnings(self) -> list[str]:
        """Warnings surfaced by profile resolution."""
        return list(self._profile.warnings)

    @property
    def haplotype_clusters(self) -> list[HaplotypeCluster]:
        """Haplotype clusters computed during the last assign() call."""
        return list(self._haplotype_clusters)

    @property
    def cn_estimate(self) -> CnEstimate | None:
        """CN estimate computed during the last assign() call."""
        return self._cn_estimate

    def assign(self, bam: BamReader, fasta: FastaReader) -> list[Assignment]:
        reads = list(build_analyzed_reads(bam, fasta, self._config, self._profile_name))

        # Haplotype clustering step
        self._haplotype_clusters = cluster_reads(reads, psv_names=self._psv_names)
        for cluster in self._haplotype_clusters:
            assign_cluster_locus(cluster, self._config)
            detect_gene_conversion(cluster, self._config)

        # Attach cluster consensus onto each read for haplotype_consistency
        read_to_cluster: dict[str, HaplotypeCluster] = {}
        for cluster in self._haplotype_clusters:
            for rid in cluster.supporting_reads:
                read_to_cluster[rid] = cluster
        for read in reads:
            maybe_cluster = read_to_cluster.get(read.read_id)
            read.cluster_consensus = (
                maybe_cluster.psv_pattern if maybe_cluster is not None else None
            )

        detected = bam.detect_tech()
        tech: ReadTech
        if detected == "ont":
            tech = "ont"
        elif detected == "short-read":
            tech = "short-read"
        else:
            tech = "ont" if reads and reads[0].is_long_read else "short-read"

        # Depth preflight pass: measure depth of primary + paralog + control regions
        depths_by_name = self._measure_depths(bam)

        # Build depths_by_locus keyed by locus IDs (for coverage_ratio adapter)
        depths_by_locus: dict[str, float] = {
            self._config.locus.id: depths_by_name[self._config.locus.id].mean_depth,
        }
        for paralog_id in self._config.coordinates.paralogs:
            paralog_stats = depths_by_name.get(paralog_id)
            depths_by_locus[paralog_id] = (
                paralog_stats.mean_depth if paralog_stats is not None else 0.0
            )

        # CN estimation
        self._cn_estimate = estimate_cn(self._config, depths_by_name, tech=tech)

        # Assemble adapters with fresh CoverageRatioEvidence
        adapters: list[EvidenceSource] = [
            *self._base_adapters,
            CoverageRatioEvidence(
                locus_id=self._config.locus.id,
                depths_by_locus=depths_by_locus,
            ),
        ]

        assignments: list[Assignment] = []
        for read in reads:
            read_cluster: HaplotypeCluster | None = read_to_cluster.get(read.read_id)
            evidences: list[EvidenceScore] = []
            for adapter in adapters:
                if not adapter.supports(tech):
                    evidences.append(
                        EvidenceScore(
                            source=adapter.name,
                            normalized=0.0,
                            raw={"reason": f"unsupported for tech={tech}"},
                            available=False,
                        )
                    )
                    continue
                evidences.append(adapter.compute([read], self._config))

            confidence, status, flags = score_assignment(
                evidences,
                self._profile,
                self._config.confidence_thresholds,
            )
            flags = set(flags)
            if read_cluster is not None and "gene_conversion_suspected" in read_cluster.notes:
                flags.add("gene_conversion_suspected")

            assignments.append(
                Assignment(
                    read_id=read.read_id,
                    assigned_locus=(
                        self._config.locus.id if status != "UNASSIGNED" else None
                    ),
                    confidence=confidence,
                    status=status,
                    evidence_scores=evidences,
                    locus_key=self._locus_key,
                    flags=flags,
                    cluster_id=read_cluster.hap_id if read_cluster is not None else None,
                )
            )
        return assignments

    def _measure_depths(self, bam: BamReader) -> dict[str, DepthStats]:
        """Measure depth for primary + paralog + control regions.

        Regions whose coordinates fall entirely outside the BAM's configured
        contig range (e.g. synthetic test BAMs that don't span full GRCh38)
        gracefully return a zero-depth DepthStats rather than propagating a
        pysam ``ValueError``.
        """
        results: dict[str, DepthStats] = {}

        primary = self._config.coordinates.primary
        results[self._config.locus.id] = self._safe_compute_depth(
            bam, primary.chrom, primary.start, primary.end,
            region_name=self._config.locus.id,
        )

        for paralog_id, coord in self._config.coordinates.paralogs.items():
            results[paralog_id] = self._safe_compute_depth(
                bam, coord.chrom, coord.start, coord.end,
                region_name=paralog_id,
            )

        for ctrl in self._config.control_regions:
            results[ctrl.name] = self._safe_compute_depth(
                bam, ctrl.chrom, ctrl.start, ctrl.end,
                region_name=ctrl.name,
            )

        return results

    @staticmethod
    def _safe_compute_depth(
        bam: BamReader,
        chrom: str,
        start: int,
        end: int,
        region_name: str,
    ) -> DepthStats:
        """Wrap ``compute_region_depth`` and fall back to zero depth when the
        region is outside the BAM's contig bounds (common for synthetic
        fixtures or contigs absent from the BAM header)."""
        try:
            return compute_region_depth(
                bam, chrom, start, end, region_name=region_name,
            )
        except (ValueError, KeyError):
            return DepthStats(
                region_name=region_name,
                chrom=chrom,
                start=start,
                end=end,
                mean_depth=0.0,
                median_depth=0.0,
                length_bp=max(0, end - start + 1),
                reads_counted=0,
            )

    def _compute_locus_key(self) -> str:
        c = self._config.coordinates.primary
        raw = f"{self._config.locus.id}:{c.chrom}:{c.start}:{c.end}"
        digest = hashlib.blake2s(raw.encode(), digest_size=3).hexdigest()
        return f"{self._config.locus.id}:{digest}"
