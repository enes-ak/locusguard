import pytest
from pydantic import ValidationError

from locusguard.config import LocusConfig

MINIMAL_VALID_CONFIG = {
    "schema_version": "1.0",
    "locusguard_compat": ">=0.1.0,<1.0.0",
    "locus": {
        "id": "SMN1",
        "name": "Survival of Motor Neuron 1",
        "gene_family": "SMN",
        "paralogs": ["SMN2"],
        "clinical_relevance": "SMA",
    },
    "reference": "grch38",
    "coordinates": {
        "primary": {"chrom": "chr5", "start": 70924941, "end": 70953012},
        "paralogs": {
            "SMN2": {"chrom": "chr5", "start": 69345350, "end": 69373421},
        },
    },
    "psvs": [
        {
            "name": "c.840C>T",
            "chrom": "chr5",
            "pos": 70951946,
            "alleles": {"SMN1": "C", "SMN2": "T"},
            "exon": 7,
            "discriminating_power": "high",
        }
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


def test_minimal_valid_config_parses():
    cfg = LocusConfig.model_validate(MINIMAL_VALID_CONFIG)
    assert cfg.locus.id == "SMN1"
    assert cfg.coordinates.primary.chrom == "chr5"
    assert len(cfg.psvs) == 1
    assert cfg.psvs[0].alleles == {"SMN1": "C", "SMN2": "T"}


def test_thresholds_require_resolved_above_probable():
    bad = {**MINIMAL_VALID_CONFIG, "confidence_thresholds": {"resolved": 0.4, "probable": 0.6}}
    with pytest.raises(ValidationError, match="resolved.*must be.*greater than.*probable"):
        LocusConfig.model_validate(bad)


def test_coord_start_must_be_less_than_end():
    bad = {**MINIMAL_VALID_CONFIG}
    bad_coords = {**bad["coordinates"]}
    bad_coords["primary"] = {"chrom": "chr5", "start": 70953012, "end": 70924941}
    bad["coordinates"] = bad_coords
    with pytest.raises(ValidationError, match="start.*less than.*end"):
        LocusConfig.model_validate(bad)


def test_missing_required_field_fails():
    bad = {k: v for k, v in MINIMAL_VALID_CONFIG.items() if k != "psvs"}
    with pytest.raises(ValidationError, match="psvs"):
        LocusConfig.model_validate(bad)


def test_locus_config_supports_control_regions():
    config_with_controls = {
        **MINIMAL_VALID_CONFIG,
        "control_regions": [
            {"name": "chr5_singlecopy_A", "chrom": "chr5", "start": 85000000, "end": 85010000},
            {"name": "chr5_singlecopy_B", "chrom": "chr5", "start": 90000000, "end": 90010000},
        ],
    }
    cfg = LocusConfig.model_validate(config_with_controls)
    assert len(cfg.control_regions) == 2
    assert cfg.control_regions[0].name == "chr5_singlecopy_A"


def test_locus_config_control_regions_optional_and_defaults_to_empty():
    cfg = LocusConfig.model_validate(MINIMAL_VALID_CONFIG)
    assert cfg.control_regions == []
