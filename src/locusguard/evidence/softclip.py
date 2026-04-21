"""Soft-clip based evidence.

Reads within their assigned locus should have minimal soft-clipping. Large
soft-clips at paralog boundaries signal misalignment — the aligner couldn't
place the tail, which may mean the tail actually belongs to the paralog copy.

Small soft-clips (<10bp total) are base-caller noise and are ignored.
"""
from __future__ import annotations

from collections.abc import Iterable

from locusguard.config.schema import LocusConfig
from locusguard.evidence.base import ReadTech
from locusguard.types import AnalyzedRead, EvidenceScore


_NOISE_THRESHOLD_BP = 10
_DEFAULT_READ_LENGTH = 2000


class SoftclipEvidence:
    name = "softclip"

    def supports(self, tech: ReadTech) -> bool:
        return True

    def compute(
        self,
        reads: Iterable[AnalyzedRead],
        locus_config: LocusConfig,
    ) -> EvidenceScore:
        reads = list(reads)
        if not reads:
            return EvidenceScore(
                source=self.name,
                normalized=0.0,
                raw={"reason": "no reads"},
                available=False,
            )

        per_read_scores: list[float] = []
        for read in reads:
            total_sc = read.softclip_5p + read.softclip_3p
            if total_sc < _NOISE_THRESHOLD_BP:
                per_read_scores.append(1.0)
                continue
            # Use default read length if unknown; real read length would come
            # from len(query_sequence), but AnalyzedRead does not carry it.
            sc_fraction = min(1.0, total_sc / _DEFAULT_READ_LENGTH)
            per_read_scores.append(1.0 - sc_fraction)

        mean = sum(per_read_scores) / len(per_read_scores)
        return EvidenceScore(
            source=self.name,
            normalized=mean,
            raw={
                "reads": len(reads),
                "mean_softclip_5p": sum(r.softclip_5p for r in reads) / len(reads),
                "mean_softclip_3p": sum(r.softclip_3p for r in reads) / len(reads),
            },
            available=True,
        )
