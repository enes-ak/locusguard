"""Haplotype consistency evidence — long-read only.

Measures whether a read's PSV observations agree with its cluster's consensus.
Requires at least 2 reach-True PSV observations for a meaningful score.
"""
from __future__ import annotations

from collections.abc import Iterable

from locusguard.config.schema import LocusConfig
from locusguard.types import AnalyzedRead, EvidenceScore

_MIN_OVERLAP = 2


class HaplotypeConsistencyEvidence:
    name = "haplotype_consistency"

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

        agree = 0
        disagree = 0
        reads_with_consensus = 0
        for read in reads:
            cons = read.cluster_consensus
            if cons is None:
                continue
            reads_with_consensus += 1
            overlap = 0
            read_agree = 0
            read_disagree = 0
            for psv_name, cons_base in cons.items():
                obs = read.psv_observations.get(psv_name)
                if obs is None or not obs.reach or cons_base == "?":
                    continue
                overlap += 1
                if obs.base == cons_base:
                    read_agree += 1
                else:
                    read_disagree += 1
            if overlap < _MIN_OVERLAP:
                # read contributes nothing — require >=2 overlapping PSVs
                continue
            agree += read_agree
            disagree += read_disagree

        total = agree + disagree
        if total == 0 or reads_with_consensus == 0:
            return EvidenceScore(
                source=self.name,
                normalized=0.0,
                raw={"reason": "insufficient cluster overlap"},
                available=False,
            )

        return EvidenceScore(
            source=self.name,
            normalized=agree / total,
            raw={"agree": agree, "disagree": disagree, "reads": reads_with_consensus},
            available=True,
        )
