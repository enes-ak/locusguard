import pytest

from locusguard.config.schema import LocusConfig
from locusguard.evidence.softclip import SoftclipEvidence
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


def _read(sc5: int, sc3: int, read_len: int = 2000) -> AnalyzedRead:
    # Fake read_length by filling psv_observations / using mapq and relying
    # on adapter to use softclip fields.
    return AnalyzedRead(
        read_id="r",
        aligned_chrom="chr5",
        aligned_pos=120,
        psv_observations={"P1": PSVObs(base="C", qual=30, reach=True)},
        mapq=60,
        softclip_5p=sc5,
        softclip_3p=sc3,
        is_long_read=True,
        is_supplementary=False,
        original_mapq_zero=False,
    )


def test_no_softclip_scores_one(cfg):
    score = SoftclipEvidence().compute([_read(0, 0)], cfg)
    assert score.available is True
    assert score.normalized == pytest.approx(1.0)


def test_noise_below_threshold_ignored(cfg):
    # 5bp each end is under the 10bp total threshold
    score = SoftclipEvidence().compute([_read(3, 6)], cfg)
    assert score.normalized == pytest.approx(1.0)


def test_heavy_softclip_scores_low(cfg):
    # 500bp soft-clip on a 2000bp read — 25% soft-clipped → score 0.75
    score = SoftclipEvidence().compute([_read(500, 0)], cfg)
    assert 0.7 < score.normalized < 0.8


def test_supports_both_techs():
    assert SoftclipEvidence().supports("ont") is True
    assert SoftclipEvidence().supports("short-read") is True


def test_empty_reads_unavailable(cfg):
    score = SoftclipEvidence().compute([], cfg)
    assert score.available is False
