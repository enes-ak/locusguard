import pytest

from locusguard.config.schema import LocusConfig
from locusguard.evidence.haplotype_consistency import HaplotypeConsistencyEvidence
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
    "psvs": [
        {"name": "P1", "chrom": "chr5", "pos": 150, "alleles": {"SMN1": "C", "SMN2": "T"}},
        {"name": "P2", "chrom": "chr5", "pos": 200, "alleles": {"SMN1": "A", "SMN2": "G"}},
    ],
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


def _read(**psvs_and_consensus) -> AnalyzedRead:
    cons = psvs_and_consensus.pop("consensus", None)
    obs = {}
    for name, base in psvs_and_consensus.items():
        obs[name] = PSVObs(
            base=base if base != "?" else "N",
            qual=30,
            reach=base != "?",
        )
    return AnalyzedRead(
        read_id="r",
        aligned_chrom="chr5",
        aligned_pos=120,
        psv_observations=obs,
        mapq=60,
        softclip_5p=0,
        softclip_3p=0,
        is_long_read=True,
        is_supplementary=False,
        original_mapq_zero=False,
        cluster_consensus=cons,
    )


def test_read_fully_consistent_with_cluster_scores_one(cfg):
    adapter = HaplotypeConsistencyEvidence()
    score = adapter.compute(
        [_read(P1="C", P2="A", consensus={"P1": "C", "P2": "A"})],
        cfg,
    )
    assert score.available is True
    assert score.normalized == pytest.approx(1.0)


def test_one_disagreement_scores_half(cfg):
    adapter = HaplotypeConsistencyEvidence()
    score = adapter.compute(
        [_read(P1="C", P2="G", consensus={"P1": "C", "P2": "A"})],
        cfg,
    )
    assert score.normalized == pytest.approx(0.5)


def test_missing_cluster_consensus_unavailable(cfg):
    adapter = HaplotypeConsistencyEvidence()
    score = adapter.compute([_read(P1="C", P2="A", consensus=None)], cfg)
    assert score.available is False


def test_unsupports_short_read():
    adapter = HaplotypeConsistencyEvidence()
    assert adapter.supports("ont") is True
    assert adapter.supports("short-read") is False


def test_only_one_reach_psv_unavailable(cfg):
    # Need >=2 overlapping observations for consistency to be meaningful.
    adapter = HaplotypeConsistencyEvidence()
    score = adapter.compute(
        [_read(P1="C", P2="?", consensus={"P1": "C", "P2": "A"})],
        cfg,
    )
    assert score.available is False
