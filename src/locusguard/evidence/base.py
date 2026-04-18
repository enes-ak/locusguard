"""Evidence source protocol.

Concrete evidence adapters implement this interface. They are the pluggable
unit of scoring: the scoring layer iterates over enabled adapters, collects
an EvidenceScore from each, and combines them per the active profile.
"""
from __future__ import annotations

from collections.abc import Iterable
from typing import Literal, Protocol, runtime_checkable

from locusguard.config.schema import LocusConfig
from locusguard.types import AnalyzedRead, EvidenceScore

ReadTech = Literal["ont", "short-read"]


@runtime_checkable
class EvidenceSource(Protocol):
    """Contract for pluggable evidence adapters."""

    name: str

    def supports(self, tech: ReadTech) -> bool:
        """Return True if this adapter produces meaningful output for `tech`."""
        ...

    def compute(
        self,
        reads: Iterable[AnalyzedRead],
        locus_config: LocusConfig,
    ) -> EvidenceScore:
        """Compute an EvidenceScore from the given reads in the given locus context.

        Implementations must return an EvidenceScore with available=False if
        the underlying signal is insufficient (e.g., too few reads observed).
        """
        ...
