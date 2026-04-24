"""MAPQ-based evidence for ONT long reads.

High MAPQ strongly supports the aligner's choice of locus. When an
individual read is not flagged as long (read.is_long_read=False, rare
in ONT pipelines), MAPQ=0 is treated as a neutral prior (0.5) rather
than a penalty — in those cases the short-read convention (BWA-MEM-style
equally-good-placement MAPQ=0) may apply, so we avoid penalizing it.
"""
from __future__ import annotations

from collections.abc import Iterable

from locusguard.config.schema import LocusConfig
from locusguard.types import AnalyzedRead, EvidenceScore

_MAX_MAPQ = 60


class MapqPatternEvidence:
    name = "mapq_pattern"

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
