"""MAPQ-based evidence — tech-differentiated.

Long-read: high MAPQ strongly supports the aligner's choice of locus.
Short-read: MAPQ=0 is common in paralog regions (BWA-MEM assigns MAPQ=0 to
equally-good placements), so it's not a strong *negative* signal. We treat
it as a neutral prior (0.5) rather than a penalty.
"""
from __future__ import annotations

from collections.abc import Iterable

from locusguard.config.schema import LocusConfig
from locusguard.evidence.base import ReadTech
from locusguard.types import AnalyzedRead, EvidenceScore

_MAX_MAPQ = 60


class MapqPatternEvidence:
    name = "mapq_pattern"

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

        scores: list[float] = []
        mapq_zero_count = 0
        for read in reads:
            if read.is_long_read:
                val = min(read.mapq / _MAX_MAPQ, 1.0)
                if read.original_mapq_zero:
                    val = 0.0
            elif read.mapq == 0:
                val = 0.5
                mapq_zero_count += 1
            else:
                val = 0.5 + 0.5 * min(read.mapq / _MAX_MAPQ, 1.0)
            scores.append(val)

        mean = sum(scores) / len(scores)
        return EvidenceScore(
            source=self.name,
            normalized=mean,
            raw={
                "reads": len(reads),
                "mapq_zero_fraction": mapq_zero_count / len(reads),
                "mean_mapq": sum(r.mapq for r in reads) / len(reads),
            },
            available=True,
        )
