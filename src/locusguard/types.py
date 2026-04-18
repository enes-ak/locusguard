"""Core data structures shared across the engine."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Literal

Status = Literal["RESOLVED", "PROBABLE", "AMBIGUOUS", "UNASSIGNED"]


@dataclass(frozen=True, slots=True)
class PSVObs:
    """A single PSV observation in a single read."""
    base: str           # Observed base (A/C/G/T/N)
    qual: int           # Base quality at the PSV position
    reach: bool         # True if the read physically covered the PSV position


@dataclass(slots=True)
class AnalyzedRead:
    """A read in the context of a locus analysis, with extracted features."""
    read_id: str
    aligned_chrom: str
    aligned_pos: int
    psv_observations: dict[str, PSVObs]
    mapq: int
    softclip_5p: int
    softclip_3p: int
    is_long_read: bool
    is_supplementary: bool
    original_mapq_zero: bool


@dataclass(frozen=True, slots=True)
class EvidenceScore:
    """Output of a single evidence adapter for a read or read group."""
    source: str                          # e.g., "psv_match"
    normalized: float                    # in [0, 1]
    raw: dict[str, Any]                  # source-specific raw values
    available: bool                      # False if adapter was disabled


@dataclass(slots=True)
class HaplotypeCluster:
    """A group of reads sharing a coherent PSV pattern."""
    hap_id: str
    supporting_reads: list[str]
    psv_pattern: dict[str, str]          # psv_name -> consensus base
    assigned_locus: str | None
    confidence: float
    notes: list[str] = field(default_factory=list)


@dataclass(slots=True)
class Assignment:
    """Final per-read locus assignment."""
    read_id: str
    assigned_locus: str | None           # None -> UNASSIGNED
    confidence: float
    status: Status
    evidence_scores: list[EvidenceScore]
    locus_key: str
    flags: set[str] = field(default_factory=set)
