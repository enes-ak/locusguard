import pytest

from locusguard.config.schema import LocusConfig
from locusguard.evidence.coverage_ratio import CoverageRatioEvidence
from locusguard.types import AnalyzedRead, PSVObs


_MINIMAL = {
    "schema_version": "1.0",
    "locusguard_compat": ">=0.1.0,<1.0.0",
    "locus": {"id": "SMN1", "name": "SMN1", "gene_family": "SMN", "paralogs": ["SMN2"]},
    "reference": "grch38",
    "coordinates": {
        "primary": {"chrom": "chr5", "start": 100, "end": 500},
        "paralogs": {"SMN2": {"chrom": "chr5", "start": 1000, "end": 1400}},
    },
    "psvs": [{"name": "P1", "chrom": "chr5", "pos": 150, "alleles": {"SMN1": "C", "SMN2": "T"}}],
    "evidence_weights": {
        "default": {
            "psv_match": 0.40, "haplotype_consistency": 0.25, "mapq_pattern": 0.10,
            "softclip": 0.05, "unique_kmer": 0.10, "coverage_ratio": 0.10,
        },
        "profile_overrides": {},
    },
    "confidence_thresholds": {"resolved": 0.80, "probable": 0.50},
}


@pytest.fixture
def cfg() -> LocusConfig:
    return LocusConfig.model_validate(_MINIMAL)


def _read() -> AnalyzedRead:
    return AnalyzedRead(
        read_id="r1",
        aligned_chrom="chr5",
        aligned_pos=120,
        psv_observations={"P1": PSVObs(base="C", qual=30, reach=True)},
        mapq=60,
        softclip_5p=0,
        softclip_3p=0,
        is_long_read=True,
        is_supplementary=False,
        original_mapq_zero=False,
    )


def test_dominant_locus_scores_high(cfg):
    adapter = CoverageRatioEvidence(
        locus_id="SMN1",
        depths_by_locus={"SMN1": 60.0, "SMN2": 30.0},
    )
    score = adapter.compute([_read()], cfg)
    assert score.source == "coverage_ratio"
    assert score.available is True
    # 60 / (60 + 30) = 0.67
    assert abs(score.normalized - (60.0 / 90.0)) < 1e-6


def test_equal_depths_scores_half(cfg):
    adapter = CoverageRatioEvidence(
        locus_id="SMN1",
        depths_by_locus={"SMN1": 30.0, "SMN2": 30.0},
    )
    score = adapter.compute([_read()], cfg)
    assert abs(score.normalized - 0.5) < 1e-6


def test_zero_total_depth_neutral(cfg):
    adapter = CoverageRatioEvidence(
        locus_id="SMN1",
        depths_by_locus={"SMN1": 0.0, "SMN2": 0.0},
    )
    score = adapter.compute([_read()], cfg)
    assert score.normalized == 0.5  # neutral
    assert score.available is True


def test_supports_ont_only():
    adapter = CoverageRatioEvidence(
        locus_id="SMN1",
        depths_by_locus={"SMN1": 30.0, "SMN2": 30.0},
    )
    assert adapter.supports("ont") is True
    assert adapter.supports("short-read") is False
