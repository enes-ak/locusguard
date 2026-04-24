"""PSV-match evidence adapter.

Given reads and a locus config, compute the fraction of PSV observations
across reads whose bases match the target locus's expected allele.
"""
from __future__ import annotations

from collections.abc import Iterable

from locusguard.config.schema import LocusConfig
from locusguard.types import AnalyzedRead, EvidenceScore


class PSVEvidence:
    """Fraction of PSV observations supporting the target locus."""

    name = "psv_match"

    def __init__(self, target_locus: str) -> None:
        self._target = target_locus

    def compute(
        self,
        reads: Iterable[AnalyzedRead],
        locus_config: LocusConfig,
    ) -> EvidenceScore:
        psv_name_to_expected = {
            p.name: p.alleles.get(self._target)
            for p in locus_config.psvs
            if self._target in p.alleles
        }

        matches = 0
        mismatches = 0
        for read in reads:
            for psv_name, expected_base in psv_name_to_expected.items():
                if expected_base is None:
                    continue
                obs = read.psv_observations.get(psv_name)
                if obs is None or not obs.reach:
                    continue
                if obs.base == expected_base:
                    matches += 1
                else:
                    mismatches += 1

        total = matches + mismatches
        if total == 0:
            return EvidenceScore(
                source=self.name,
                normalized=0.0,
                raw={"matches": 0, "mismatches": 0, "total_observations": 0},
                available=False,
            )

        normalized = matches / total
        return EvidenceScore(
            source=self.name,
            normalized=normalized,
            raw={
                "matches": matches,
                "mismatches": mismatches,
                "total_observations": total,
                "target_locus": self._target,
            },
            available=True,
        )
