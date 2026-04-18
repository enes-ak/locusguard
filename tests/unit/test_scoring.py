import pytest

from locusguard.config.resolver import ResolvedProfile
from locusguard.config.schema import EvidenceWeights, Thresholds
from locusguard.scoring import score_assignment
from locusguard.types import EvidenceScore


def _weights(**kw) -> EvidenceWeights:
    base = {
        "psv_match": 0.40,
        "haplotype_consistency": 0.25,
        "mapq_pattern": 0.10,
        "softclip": 0.05,
        "unique_kmer": 0.10,
        "coverage_ratio": 0.10,
    }
    base.update(kw)
    return EvidenceWeights(**base)


def _profile(**kw) -> ResolvedProfile:
    weights = kw.pop("weights", _weights())
    cap = kw.pop("cap_confidence", 1.0)
    disabled = kw.pop("disabled", set())
    return ResolvedProfile(
        weights=weights,
        cap_confidence=cap,
        disabled=disabled,
    )


_DEFAULT_THRESHOLDS = Thresholds(resolved=0.80, probable=0.50)


def test_single_high_evidence_yields_resolved():
    evidences = [EvidenceScore(source="psv_match", normalized=0.95, raw={}, available=True)]
    conf, status, flags = score_assignment(evidences, _profile(), _DEFAULT_THRESHOLDS)
    assert conf == pytest.approx(0.95)
    assert status == "RESOLVED"
    assert flags == set()


def test_weighted_average_combines_multiple_evidences():
    evidences = [
        EvidenceScore(source="psv_match", normalized=1.0, raw={}, available=True),
        EvidenceScore(source="haplotype_consistency", normalized=0.0, raw={}, available=True),
    ]
    # weights: psv_match=0.40, haplotype_consistency=0.25 -> denom = 0.65
    # score = (1.0*0.40 + 0.0*0.25) / 0.65 = 0.615
    conf, status, _ = score_assignment(evidences, _profile(), _DEFAULT_THRESHOLDS)
    assert conf == pytest.approx(0.40 / 0.65, abs=1e-6)
    assert status == "PROBABLE"


def test_unavailable_evidence_excluded_from_average():
    evidences = [
        EvidenceScore(source="psv_match", normalized=0.9, raw={}, available=True),
        EvidenceScore(source="unique_kmer", normalized=0.0, raw={}, available=False),
    ]
    conf, status, _ = score_assignment(evidences, _profile(), _DEFAULT_THRESHOLDS)
    # Only psv_match counts -> 0.9
    assert conf == pytest.approx(0.9)
    assert status == "RESOLVED"


def test_cap_confidence_clips_score():
    evidences = [EvidenceScore(source="psv_match", normalized=1.0, raw={}, available=True)]
    profile = _profile(cap_confidence=0.85)
    conf, status, _ = score_assignment(evidences, profile, _DEFAULT_THRESHOLDS)
    assert conf == pytest.approx(0.85)
    assert status == "RESOLVED"


def test_disabled_evidence_excluded():
    evidences = [
        EvidenceScore(source="psv_match", normalized=0.9, raw={}, available=True),
        EvidenceScore(source="coverage_ratio", normalized=0.2, raw={}, available=True),
    ]
    profile = _profile(disabled={"coverage_ratio"})
    conf, _, _ = score_assignment(evidences, profile, _DEFAULT_THRESHOLDS)
    # Only psv_match contributes
    assert conf == pytest.approx(0.9)


def test_no_available_evidence_yields_unassigned():
    evidences = [EvidenceScore(source="psv_match", normalized=0.0, raw={}, available=False)]
    conf, status, flags = score_assignment(evidences, _profile(), _DEFAULT_THRESHOLDS)
    assert conf == 0.0
    assert status == "UNASSIGNED"
    assert "insufficient_coverage" in flags


def test_low_score_yields_ambiguous():
    evidences = [EvidenceScore(source="psv_match", normalized=0.3, raw={}, available=True)]
    conf, status, _ = score_assignment(evidences, _profile(), _DEFAULT_THRESHOLDS)
    assert status == "AMBIGUOUS"
