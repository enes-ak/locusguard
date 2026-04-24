import pytest

from locusguard.config.schema import LocusConfig
from locusguard.evidence.mapq_pattern import MapqPatternEvidence
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


def _read(mapq: int, mapq_zero: bool = False) -> AnalyzedRead:
    return AnalyzedRead(
        read_id="r",
        aligned_chrom="chr5",
        aligned_pos=120,
        psv_observations={"P1": PSVObs(base="C", qual=30, reach=True)},
        mapq=mapq,
        softclip_5p=0,
        softclip_3p=0,
        is_long_read=True,
        is_supplementary=False,
        original_mapq_zero=mapq_zero,
    )


def test_ont_high_mapq_scores_high(cfg):
    adapter = MapqPatternEvidence()
    score = adapter.compute([_read(60)], cfg)
    assert score.source == "mapq_pattern"
    assert score.available is True
    assert score.normalized >= 0.99


def test_ont_zero_mapq_scores_zero(cfg):
    adapter = MapqPatternEvidence()
    score = adapter.compute([_read(0, mapq_zero=True)], cfg)
    assert score.normalized == 0.0


def test_empty_reads_unavailable(cfg):
    adapter = MapqPatternEvidence()
    score = adapter.compute([], cfg)
    assert score.available is False
