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

from locusguard.config.resolver import ResolvedProfile, resolve_profile
from locusguard.config.schema import PSV, LocusConfig
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
        self._adapters = [
            PSVEvidence(target_locus=config.locus.id),
            MapqPatternEvidence(),
            SoftclipEvidence(),
            HaplotypeConsistencyEvidence(),
        ]
        self._haplotype_clusters: list[HaplotypeCluster] = []

    @property
    def warnings(self) -> list[str]:
        """Warnings surfaced by profile resolution."""
        return list(self._profile.warnings)

    @property
    def haplotype_clusters(self) -> list[HaplotypeCluster]:
        """Haplotype clusters computed during the last assign() call."""
        return list(self._haplotype_clusters)

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
            cluster = read_to_cluster.get(read.read_id)
            read.cluster_consensus = cluster.psv_pattern if cluster is not None else None

        tech = "ont" if reads and reads[0].is_long_read else "short-read"
        assignments: list[Assignment] = []
        for read in reads:
            cluster = read_to_cluster.get(read.read_id)
            evidences: list[EvidenceScore] = []
            for adapter in self._adapters:
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
            if cluster is not None and "gene_conversion_suspected" in cluster.notes:
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
                    cluster_id=cluster.hap_id if cluster is not None else None,
                )
            )
        return assignments

    def _compute_locus_key(self) -> str:
        c = self._config.coordinates.primary
        raw = f"{self._config.locus.id}:{c.chrom}:{c.start}:{c.end}"
        digest = hashlib.blake2s(raw.encode(), digest_size=3).hexdigest()
        return f"{self._config.locus.id}:{digest}"
