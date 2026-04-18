"""Pydantic models for LocusGuard locus configuration YAML."""
from __future__ import annotations

from typing import Any, Literal

from pydantic import BaseModel, ConfigDict, Field, model_validator


class CoordRange(BaseModel):
    """A genomic coordinate range."""
    model_config = ConfigDict(extra="forbid")

    chrom: str
    start: int = Field(ge=0)
    end: int = Field(ge=1)

    @model_validator(mode="after")
    def _check_order(self) -> CoordRange:
        if self.start >= self.end:
            raise ValueError("start must be less than end")
        return self


class Coordinates(BaseModel):
    model_config = ConfigDict(extra="forbid")

    primary: CoordRange
    paralogs: dict[str, CoordRange] = Field(default_factory=dict)


class PSV(BaseModel):
    """Paralog-Specific Variant — a position that distinguishes paralogs."""
    model_config = ConfigDict(extra="forbid")

    name: str
    chrom: str
    pos: int = Field(ge=0)
    alleles: dict[str, str]
    exon: int | None = None
    discriminating_power: Literal["low", "medium", "high"] = "high"
    notes: str | None = None


class LocusMetadata(BaseModel):
    model_config = ConfigDict(extra="forbid")

    id: str
    name: str
    gene_family: str
    paralogs: list[str]
    clinical_relevance: str | None = None


class CriticalRegion(BaseModel):
    model_config = ConfigDict(extra="forbid")

    name: str
    chrom: str
    start: int
    end: int
    required_coverage: int = 10


class GeneConversion(BaseModel):
    model_config = ConfigDict(extra="forbid")

    known_hotspots: list[CriticalRegion] = Field(default_factory=list)


class CopyNumber(BaseModel):
    model_config = ConfigDict(extra="forbid")

    expected_total: int | None = None
    individual_range: tuple[int, int] | None = None
    estimation_method: str = "depth_ratio"


class EvidenceWeights(BaseModel):
    """Default weights; profile overrides applied at resolution time."""
    model_config = ConfigDict(extra="forbid")

    psv_match: float = Field(ge=0)
    haplotype_consistency: float = Field(ge=0)
    mapq_pattern: float = Field(ge=0)
    softclip: float = Field(ge=0)
    unique_kmer: float = Field(ge=0)
    coverage_ratio: float = Field(ge=0)


class ProfileOverride(BaseModel):
    """Per-profile override on top of default weights."""
    model_config = ConfigDict(extra="allow")  # allow weight field overrides dynamically

    cap_confidence: float = Field(default=1.0, ge=0, le=1)
    enable: list[str] = Field(default_factory=list)
    disable: list[str] = Field(default_factory=list)


class Thresholds(BaseModel):
    model_config = ConfigDict(extra="forbid")

    resolved: float = Field(ge=0, le=1)
    probable: float = Field(ge=0, le=1)

    @model_validator(mode="after")
    def _check_order(self) -> Thresholds:
        if self.resolved <= self.probable:
            raise ValueError("resolved threshold must be greater than probable threshold")
        return self


class ProfiledWeights(BaseModel):
    model_config = ConfigDict(extra="forbid")

    default: EvidenceWeights
    profile_overrides: dict[str, ProfileOverride] = Field(default_factory=dict)


class ValidationSample(BaseModel):
    model_config = ConfigDict(extra="forbid")

    sample: str
    expected: dict[str, int | str]
    source: str | None = None


class ValidationBlock(BaseModel):
    model_config = ConfigDict(extra="forbid")

    tested_samples: list[ValidationSample] = Field(default_factory=list)
    tested_tech: list[str] = Field(default_factory=list)
    config_author: str | None = None
    last_calibrated: str | None = None


class ReferenceCitation(BaseModel):
    model_config = ConfigDict(extra="forbid")

    doi: str | None = None
    note: str | None = None


class WarningRule(BaseModel):
    model_config = ConfigDict(extra="forbid")

    profile: str | None = None
    condition: str | None = None
    always: str | None = None
    message: str | None = None


class LocusConfig(BaseModel):
    """Top-level locus configuration."""
    model_config = ConfigDict(extra="forbid")

    schema_version: str
    locusguard_compat: str
    locus: LocusMetadata
    reference: str
    coordinates: Coordinates
    critical_regions: list[CriticalRegion] = Field(default_factory=list)
    psvs: list[PSV]
    unique_kmers: dict[str, Any] | None = None
    gene_conversion: GeneConversion | None = None
    copy_number: CopyNumber | None = None
    evidence_weights: ProfiledWeights
    confidence_thresholds: Thresholds
    warnings: list[WarningRule] = Field(default_factory=list)
    validation: ValidationBlock | None = None
    references: list[ReferenceCitation] = Field(default_factory=list)
