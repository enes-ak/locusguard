"""Per-locus assignment orchestrator.

For each locus, fetch reads overlapping the primary + paralog regions,
build AnalyzedRead objects (with PSV base extraction), run evidence
adapters, score, and emit Assignments.
"""
from __future__ import annotations

import hashlib
from collections.abc import Iterator

import pysam

from locusguard.config.resolver import ResolvedProfile, resolve_profile
from locusguard.config.schema import PSV, LocusConfig
from locusguard.evidence.psv import PSVEvidence
from locusguard.io.bam import BamReader
from locusguard.io.fasta import FastaReader
from locusguard.scoring import score_assignment
from locusguard.types import AnalyzedRead, Assignment, EvidenceScore, PSVObs

# pysam CIGAR op codes (BAM spec): operation integer for soft-clip.
_CIGAR_SOFT_CLIP = 4


def build_analyzed_reads(
    bam: BamReader,
    fasta: FastaReader,
    config: LocusConfig,
    profile_name: str | None,
) -> Iterator[AnalyzedRead]:
    """Iterate reads overlapping primary + paralog regions, build AnalyzedReads."""
    regions = [config.coordinates.primary, *config.coordinates.paralogs.values()]
    long_read_hint = bam.estimated_is_long_read()
    seen: set[str] = set()

    for region in regions:
        for read in bam.fetch(region.chrom, region.start, region.end):
            if read.query_name is None or read.query_name in seen:
                continue
            seen.add(read.query_name)
            analyzed = _analyze_read(read, config.psvs, long_read_hint)
            # Only yield reads that physically reach at least one configured PSV;
            # reads that reach none of the PSVs provide no information for the
            # PSV-match adapter and would be flagged UNASSIGNED with no signal.
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
        # Config PSV positions use the 1-based genomics convention (matches VCF
        # POS, literature, samtools). pysam's get_aligned_pairs returns 0-based
        # reference coordinates, so convert.
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
    )


def _ref_to_query_index(read: pysam.AlignedSegment) -> dict[int, int]:
    pairs = read.get_aligned_pairs(matches_only=True)  # [(qpos, refpos), ...]
    return {ref: q for q, ref in pairs}


def _softclip_amounts(read: pysam.AlignedSegment) -> tuple[int, int]:
    cigar = read.cigartuples or []
    if not cigar:
        return 0, 0
    sc_5 = cigar[0][1] if cigar[0][0] == _CIGAR_SOFT_CLIP else 0
    sc_3 = cigar[-1][1] if cigar[-1][0] == _CIGAR_SOFT_CLIP else 0
    return sc_5, sc_3


class LocusAssigner:
    """Orchestrates evidence collection and scoring for a single locus."""

    def __init__(self, config: LocusConfig, profile_name: str | None) -> None:
        self._config = config
        self._profile_name = profile_name
        self._profile: ResolvedProfile = resolve_profile(
            config.evidence_weights, profile_name
        )
        # Phase 1: only the PSV evidence adapter is wired up
        self._adapters = [PSVEvidence(target_locus=config.locus.id)]
        self._locus_key = self._compute_locus_key()

    @property
    def warnings(self) -> list[str]:
        """Warnings surfaced by profile resolution (e.g., unknown profile name)."""
        return list(self._profile.warnings)

    def assign(self, bam: BamReader, fasta: FastaReader) -> list[Assignment]:
        assignments: list[Assignment] = []
        for read in build_analyzed_reads(bam, fasta, self._config, self._profile_name):
            evidences: list[EvidenceScore] = []
            for adapter in self._adapters:
                evidences.append(adapter.compute([read], self._config))
            confidence, status, flags = score_assignment(
                evidences,
                self._profile,
                self._config.confidence_thresholds,
            )
            assignments.append(
                Assignment(
                    read_id=read.read_id,
                    assigned_locus=self._config.locus.id if status != "UNASSIGNED" else None,
                    confidence=confidence,
                    status=status,
                    evidence_scores=evidences,
                    locus_key=self._locus_key,
                    flags=flags,
                )
            )
        return assignments

    def _compute_locus_key(self) -> str:
        c = self._config.coordinates.primary
        raw = f"{self._config.locus.id}:{c.chrom}:{c.start}:{c.end}"
        digest = hashlib.blake2s(raw.encode(), digest_size=3).hexdigest()
        return f"{self._config.locus.id}:{digest}"
