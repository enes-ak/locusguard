import pytest

from locusguard.config.schema import LocusConfig
from locusguard.haplotype.consensus import (
    assign_cluster_locus,
    detect_gene_conversion,
)
from locusguard.types import HaplotypeCluster


_CFG_DICT = {
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
        {"name": "P3", "chrom": "chr5", "pos": 250, "alleles": {"SMN1": "T", "SMN2": "A"}},
    ],
    "evidence_weights": {
        "default": {
            "psv_match": 0.40, "haplotype_consistency": 0.25, "mapq_pattern": 0.10,
            "softclip": 0.05, "unique_kmer": 0.10, "coverage_ratio": 0.10,
        },
        "profile_overrides": {},
    },
    "confidence_thresholds": {"resolved": 0.80, "probable": 0.50},
    "gene_conversion": {
        "known_hotspots": [
            {"name": "test_hotspot", "chrom": "chr5", "start": 140, "end": 210},
        ],
    },
}


@pytest.fixture
def cfg() -> LocusConfig:
    return LocusConfig.model_validate(_CFG_DICT)


def _cluster(pattern: dict[str, str]) -> HaplotypeCluster:
    return HaplotypeCluster(
        hap_id="H1",
        supporting_reads=["r1", "r2"],
        psv_pattern=pattern,
        assigned_locus=None,
        confidence=0.0,
        notes=[],
    )


def test_pure_smn1_pattern_assigned_smn1(cfg):
    cluster = _cluster({"P1": "C", "P2": "A", "P3": "T"})
    assign_cluster_locus(cluster, cfg)
    assert cluster.assigned_locus == "SMN1"
    assert cluster.confidence >= 0.99


def test_pure_smn2_pattern_assigned_smn2(cfg):
    cluster = _cluster({"P1": "T", "P2": "G", "P3": "A"})
    assign_cluster_locus(cluster, cfg)
    assert cluster.assigned_locus == "SMN2"


def test_mixed_pattern_flags_gene_conversion(cfg):
    cluster = _cluster({"P1": "C", "P2": "A", "P3": "A"})  # P3 is SMN2-like
    assign_cluster_locus(cluster, cfg)
    detect_gene_conversion(cluster, cfg)
    assert "gene_conversion_suspected" in cluster.notes


def test_mixed_pattern_in_known_hotspot_gets_hotspot_note(cfg):
    # P1 (pos 150) and P2 (pos 200) are within hotspot 140-210.
    # Construct a cluster where P1 is SMN1-like and P2 is SMN2-like — mixed inside hotspot.
    cluster = _cluster({"P1": "C", "P2": "G", "P3": "T"})
    assign_cluster_locus(cluster, cfg)
    detect_gene_conversion(cluster, cfg)
    assert "gene_conversion_suspected" in cluster.notes
    assert "hotspot_match:test_hotspot" in cluster.notes


def test_all_unknown_pattern_unassigned(cfg):
    cluster = _cluster({"P1": "?", "P2": "?", "P3": "?"})
    assign_cluster_locus(cluster, cfg)
    assert cluster.assigned_locus is None


def test_partial_observation_still_assignable(cfg):
    cluster = _cluster({"P1": "C", "P2": "?", "P3": "T"})
    assign_cluster_locus(cluster, cfg)
    assert cluster.assigned_locus == "SMN1"
    assert cluster.confidence < 1.0  # missing P2 reduces confidence
