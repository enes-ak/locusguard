import pytest

from locusguard.config.schema import LocusConfig
from locusguard.evidence.psv import PSVEvidence
from locusguard.types import AnalyzedRead, PSVObs

MINIMAL_SMN1_CFG_DICT = {
    "schema_version": "1.0",
    "locusguard_compat": ">=0.1.0,<1.0.0",
    "locus": {
        "id": "SMN1",
        "name": "SMN1",
        "gene_family": "SMN",
        "paralogs": ["SMN2"],
    },
    "reference": "grch38",
    "coordinates": {
        "primary": {"chrom": "chr5", "start": 12000, "end": 15000},
        "paralogs": {
            "SMN2": {"chrom": "chr5", "start": 2000, "end": 5000},
        },
    },
    "psvs": [
        {
            "name": "c.840C>T",
            "chrom": "chr5",
            "pos": 14000,
            "alleles": {"SMN1": "C", "SMN2": "T"},
            "exon": 7,
        },
        {
            "name": "p2",
            "chrom": "chr5",
            "pos": 14500,
            "alleles": {"SMN1": "A", "SMN2": "G"},
        },
    ],
    "evidence_weights": {
        "default": {
            "psv_match": 0.40,
            "haplotype_consistency": 0.25,
            "mapq_pattern": 0.10,
            "softclip": 0.05,
            "unique_kmer": 0.10,
            "coverage_ratio": 0.10,
        },
        "profile_overrides": {},
    },
    "confidence_thresholds": {"resolved": 0.80, "probable": 0.50},
}


@pytest.fixture
def smn1_cfg() -> LocusConfig:
    return LocusConfig.model_validate(MINIMAL_SMN1_CFG_DICT)


def _read_with_psvs(**psvs) -> AnalyzedRead:
    return AnalyzedRead(
        read_id="r",
        aligned_chrom="chr5",
        aligned_pos=13000,
        psv_observations={
            name: PSVObs(base=base, qual=30, reach=True)
            for name, base in psvs.items()
        },
        mapq=60,
        softclip_5p=0,
        softclip_3p=0,
        is_long_read=True,
        is_supplementary=False,
        original_mapq_zero=False,
    )


def test_psv_evidence_full_match_scores_high(smn1_cfg):
    adapter = PSVEvidence(target_locus="SMN1")
    reads = [
        _read_with_psvs(**{"c.840C>T": "C", "p2": "A"}),
        _read_with_psvs(**{"c.840C>T": "C", "p2": "A"}),
    ]
    score = adapter.compute(reads, smn1_cfg)
    assert score.source == "psv_match"
    assert score.available is True
    assert score.normalized >= 0.9      # all reads agree on SMN1


def test_psv_evidence_full_mismatch_scores_zero(smn1_cfg):
    adapter = PSVEvidence(target_locus="SMN1")
    reads = [
        _read_with_psvs(**{"c.840C>T": "T", "p2": "G"}),  # SMN2-consistent
        _read_with_psvs(**{"c.840C>T": "T", "p2": "G"}),
    ]
    score = adapter.compute(reads, smn1_cfg)
    assert score.normalized <= 0.1


def test_psv_evidence_no_reach_reduces_availability(smn1_cfg):
    adapter = PSVEvidence(target_locus="SMN1")
    reads = [
        AnalyzedRead(
            read_id="r",
            aligned_chrom="chr5",
            aligned_pos=13000,
            psv_observations={
                "c.840C>T": PSVObs(base="N", qual=0, reach=False),
                "p2": PSVObs(base="N", qual=0, reach=False),
            },
            mapq=60,
            softclip_5p=0,
            softclip_3p=0,
            is_long_read=True,
            is_supplementary=False,
            original_mapq_zero=False,
        )
    ]
    score = adapter.compute(reads, smn1_cfg)
    assert score.available is False


def test_psv_evidence_reports_decomposition(smn1_cfg):
    adapter = PSVEvidence(target_locus="SMN1")
    reads = [_read_with_psvs(**{"c.840C>T": "C", "p2": "A"})]
    score = adapter.compute(reads, smn1_cfg)
    assert "matches" in score.raw
    assert "total_observations" in score.raw
