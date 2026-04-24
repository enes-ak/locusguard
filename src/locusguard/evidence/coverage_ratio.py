"""Per-locus uniform coverage ratio evidence.

Score = depth(assigned_locus) / (sum depth across all configured paralogs).
Uniform across all reads from the same locus; depth is a region-level
property, not a per-read property. Complements per-read signals (PSV, MAPQ,
softclip, haplotype_consistency) by answering: "does the assigned locus
actually have depth to spare for this read?"
"""
from __future__ import annotations

from collections.abc import Iterable

from locusguard.config.schema import LocusConfig
from locusguard.types import AnalyzedRead, EvidenceScore


class CoverageRatioEvidence:
    name = "coverage_ratio"

    def __init__(
        self,
        locus_id: str,
        depths_by_locus: dict[str, float],
    ) -> None:
        self._locus_id = locus_id
        self._depths = depths_by_locus

    def compute(
        self,
        reads: Iterable[AnalyzedRead],
        locus_config: LocusConfig,
    ) -> EvidenceScore:
        reads_list = list(reads)
        total = sum(self._depths.values())
        locus_depth = self._depths.get(self._locus_id, 0.0)
        if total <= 0:
            return EvidenceScore(
                source=self.name,
                normalized=0.5,
                raw={
                    "locus_depth": 0.0,
                    "total_depth": 0.0,
                    "reason": "no_depth",
                    "reads": len(reads_list),
                },
                available=True,
            )
        ratio = locus_depth / total
        return EvidenceScore(
            source=self.name,
            normalized=ratio,
            raw={
                "locus_depth": locus_depth,
                "total_depth": total,
                "reads": len(reads_list),
            },
            available=True,
        )
