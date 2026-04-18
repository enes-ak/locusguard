"""Weighted-sum scoring and tier mapping for evidence lists."""
from __future__ import annotations

from locusguard.config.resolver import ResolvedProfile
from locusguard.config.schema import Thresholds
from locusguard.types import EvidenceScore, Status


def score_assignment(
    evidences: list[EvidenceScore],
    profile: ResolvedProfile,
    thresholds: Thresholds,
) -> tuple[float, Status, set[str]]:
    """Combine evidence scores -> (confidence, status, flags).

    - Skips unavailable evidences (adapter returned available=False)
    - Skips evidences whose source is in profile.disabled
    - Uses profile.weights for the weighted average
    - Clips the final score to profile.cap_confidence
    - Returns UNASSIGNED + flag if no usable evidence
    """
    weight_by_source = _weights_as_dict(profile)
    usable: list[tuple[float, float]] = []
    for ev in evidences:
        if not ev.available:
            continue
        if ev.source in profile.disabled:
            continue
        w = weight_by_source.get(ev.source, 0.0)
        if w <= 0:
            continue
        usable.append((ev.normalized, w))

    if not usable:
        return 0.0, "UNASSIGNED", {"insufficient_coverage"}

    numerator = sum(v * w for v, w in usable)
    denominator = sum(w for _, w in usable)
    score = numerator / denominator if denominator > 0 else 0.0
    score = min(score, profile.cap_confidence)

    status = _tier_for(score, thresholds)
    return score, status, set()


def _weights_as_dict(profile: ResolvedProfile) -> dict[str, float]:
    w = profile.weights
    return {
        "psv_match": w.psv_match,
        "haplotype_consistency": w.haplotype_consistency,
        "mapq_pattern": w.mapq_pattern,
        "softclip": w.softclip,
        "unique_kmer": w.unique_kmer,
        "coverage_ratio": w.coverage_ratio,
    }


def _tier_for(score: float, thresholds: Thresholds) -> Status:
    if score >= thresholds.resolved:
        return "RESOLVED"
    if score >= thresholds.probable:
        return "PROBABLE"
    return "AMBIGUOUS"
