"""Resolve locus config evidence weights for a specific runtime profile."""
from __future__ import annotations

from dataclasses import dataclass, field

from locusguard.config.schema import (
    EvidenceWeights,
    ProfiledWeights,
)


@dataclass(slots=True)
class ResolvedProfile:
    """The concrete weights + flags produced by applying a profile override."""
    weights: EvidenceWeights
    cap_confidence: float
    enabled: set[str] = field(default_factory=set)
    disabled: set[str] = field(default_factory=set)
    warnings: list[str] = field(default_factory=list)


_WEIGHT_FIELDS = {
    "psv_match",
    "haplotype_consistency",
    "mapq_pattern",
    "softclip",
    "unique_kmer",
    "coverage_ratio",
}


def resolve_profile(
    profiled: ProfiledWeights,
    profile_name: str | None,
) -> ResolvedProfile:
    """Apply a profile override on top of defaults.

    If `profile_name` is None or not defined, returns the default with cap_confidence=1.0.
    Unknown profile names emit a warning and fall back to default.
    """
    warnings: list[str] = []
    if profile_name is None:
        return ResolvedProfile(weights=profiled.default, cap_confidence=1.0)

    override = profiled.profile_overrides.get(profile_name)
    if override is None:
        warnings.append(f"profile '{profile_name}' not defined; falling back to default")
        return ResolvedProfile(
            weights=profiled.default,
            cap_confidence=1.0,
            warnings=warnings,
        )

    # Start with default weights, apply field overrides.
    default_dict = profiled.default.model_dump()
    override_extras = override.model_dump()
    for key, value in override_extras.items():
        if key in _WEIGHT_FIELDS:
            default_dict[key] = float(value)

    return ResolvedProfile(
        weights=EvidenceWeights(**default_dict),
        cap_confidence=override.cap_confidence,
        enabled=set(override.enable),
        disabled=set(override.disable),
        warnings=warnings,
    )
