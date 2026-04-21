import dataclasses

import pytest

from locusguard.types import (
    AnalyzedRead,
    Assignment,
    EvidenceScore,
    HaplotypeCluster,
    PSVObs,
)


def test_psv_obs_holds_base_and_qual():
    obs = PSVObs(base="C", qual=30, reach=True)
    assert obs.base == "C"
    assert obs.qual == 30
    assert obs.reach is True


def test_psv_obs_is_frozen():
    obs = PSVObs(base="C", qual=30, reach=True)
    with pytest.raises(dataclasses.FrozenInstanceError):
        obs.base = "T"  # type: ignore[misc]


def test_analyzed_read_holds_read_context():
    read = AnalyzedRead(
        read_id="read_001",
        aligned_chrom="chr5",
        aligned_pos=70951946,
        psv_observations={
            "c.840C>T": PSVObs(base="C", qual=30, reach=True),
        },
        mapq=60,
        softclip_5p=0,
        softclip_3p=5,
        is_long_read=True,
        is_supplementary=False,
        original_mapq_zero=False,
    )
    assert read.read_id == "read_001"
    assert read.psv_observations["c.840C>T"].base == "C"
    assert read.is_long_read is True


def test_evidence_score_is_frozen_and_has_metadata():
    es = EvidenceScore(
        source="psv_match",
        normalized=0.75,
        raw={"matches": 3, "total": 4},
        available=True,
    )
    assert es.source == "psv_match"
    assert 0.0 <= es.normalized <= 1.0
    assert es.available is True
    with pytest.raises(dataclasses.FrozenInstanceError):
        es.normalized = 0.0  # type: ignore[misc]


def test_assignment_captures_decision_and_flags():
    a = Assignment(
        read_id="read_001",
        assigned_locus="SMN1",
        confidence=0.87,
        status="RESOLVED",
        evidence_scores=[
            EvidenceScore(source="psv_match", normalized=0.87, raw={}, available=True),
        ],
        locus_key="SMN1:a3f9c2",
        flags=set(),
    )
    assert a.assigned_locus == "SMN1"
    assert a.status == "RESOLVED"


def test_haplotype_cluster_supports_note_accumulation():
    hc = HaplotypeCluster(
        hap_id="H1",
        supporting_reads=["r1", "r2"],
        psv_pattern={"c.840C>T": "C"},
        assigned_locus="SMN1",
        confidence=0.9,
        notes=["high_coverage"],
    )
    hc.notes.append("gene_conversion_suspected")
    assert "gene_conversion_suspected" in hc.notes


def test_assignment_supports_optional_cluster_id():
    a = Assignment(
        read_id="r1",
        assigned_locus="SMN1",
        confidence=0.9,
        status="RESOLVED",
        evidence_scores=[],
        locus_key="SMN1:abc",
        flags=set(),
        cluster_id="H1",
    )
    assert a.cluster_id == "H1"


def test_assignment_cluster_id_defaults_to_none():
    a = Assignment(
        read_id="r1",
        assigned_locus="SMN1",
        confidence=0.9,
        status="RESOLVED",
        evidence_scores=[],
        locus_key="SMN1:abc",
        flags=set(),
    )
    assert a.cluster_id is None
