"""Evidence source protocol.

Concrete evidence adapters implement this interface. The scoring layer
iterates over adapters, collects an EvidenceScore from each, and combines
them per the active profile.
"""
from __future__ import annotations

from collections.abc import Iterable
from typing import Protocol, runtime_checkable

from locusguard.config.schema import LocusConfig
from locusguard.types import AnalyzedRead, EvidenceScore


@runtime_checkable
class EvidenceSource(Protocol):
    """Contract for pluggable evidence adapters."""

    name: str

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
