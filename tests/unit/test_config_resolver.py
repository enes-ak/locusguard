from locusguard.config import (
    resolve_profile,
)
from locusguard.config.schema import ProfiledWeights


def _make_weights() -> ProfiledWeights:
    return ProfiledWeights.model_validate({
        "default": {
            "psv_match": 0.40,
            "haplotype_consistency": 0.25,
            "mapq_pattern": 0.10,
            "softclip": 0.05,
            "unique_kmer": 0.10,
            "coverage_ratio": 0.10,
        },
        "profile_overrides": {
            "ont_wgs": {
                "psv_match": 0.35,
                "haplotype_consistency": 0.30,
                "coverage_ratio": 0.10,
                "cap_confidence": 1.0,
                "enable": ["long_haplotype"],
            },
            "short_read_wes": {
                "psv_match": 0.50,
                "cap_confidence": 0.85,
                "disable": ["coverage_ratio", "long_haplotype"],
            },
        },
    })


def test_resolve_profile_returns_default_when_no_override():
    pw = _make_weights()
    resolved = resolve_profile(pw, profile_name=None)
    assert resolved.weights.psv_match == 0.40
    assert resolved.cap_confidence == 1.0
    assert resolved.disabled == set()


def test_resolve_profile_applies_overrides():
    pw = _make_weights()
    resolved = resolve_profile(pw, profile_name="ont_wgs")
    assert resolved.weights.psv_match == 0.35            # overridden
    assert resolved.weights.haplotype_consistency == 0.30
    assert resolved.weights.softclip == 0.05             # not in override, from default
    assert resolved.cap_confidence == 1.0
    assert "long_haplotype" in resolved.enabled


def test_resolve_profile_disables_evidences():
    pw = _make_weights()
    resolved = resolve_profile(pw, profile_name="short_read_wes")
    assert resolved.cap_confidence == 0.85
    assert "coverage_ratio" in resolved.disabled
    assert "long_haplotype" in resolved.disabled


def test_resolve_profile_unknown_profile_falls_back_to_default():
    pw = _make_weights()
    resolved = resolve_profile(pw, profile_name="mars_rover")
    assert resolved.weights.psv_match == 0.40
    assert resolved.warnings == ["profile 'mars_rover' not defined; falling back to default"]
