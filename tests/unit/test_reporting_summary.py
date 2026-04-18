import json

from locusguard.reporting.summary import write_summary
from locusguard.types import Assignment, EvidenceScore


def _make_assignment(locus: str, status: str, confidence: float) -> Assignment:
    return Assignment(
        read_id="r",
        assigned_locus=locus,
        confidence=confidence,
        status=status,
        evidence_scores=[
            EvidenceScore(
                source="psv_match",
                normalized=confidence,
                raw={},
                available=True,
            )
        ],
        locus_key=f"{locus}:a3f9c2",
        flags=set(),
    )


def test_summary_json_has_expected_top_level_fields(tmp_path):
    out = tmp_path / "summary.json"
    write_summary(
        output_path=out,
        sample_name="HG002",
        reference="grch38",
        tech="ont",
        data_type="wgs",
        runtime_seconds=42.5,
        assignments_by_locus={
            "SMN1": [_make_assignment("SMN1", "RESOLVED", 0.9)],
            "SMN2": [_make_assignment("SMN2", "RESOLVED", 0.85)],
        },
        variant_counts_by_locus={"SMN1": 3, "SMN2": 2},
    )
    doc = json.loads(out.read_text())
    assert doc["sample"] == "HG002"
    assert doc["tech"] == "ont"
    assert doc["data_type"] == "wgs"
    assert doc["reference"] == "grch38"
    assert doc["runtime_sec"] == 42.5
    assert "loci" in doc
    assert doc["loci"]["SMN1"]["status"] == "RESOLVED"
    assert doc["loci"]["SMN1"]["variants_annotated"] == 3


def test_summary_reports_mixed_status_as_dominant(tmp_path):
    out = tmp_path / "summary.json"
    write_summary(
        output_path=out,
        sample_name="HG002",
        reference="grch38",
        tech="ont",
        data_type="wgs",
        runtime_seconds=1.0,
        assignments_by_locus={
            "SMN1": [
                _make_assignment("SMN1", "RESOLVED", 0.9),
                _make_assignment("SMN1", "RESOLVED", 0.9),
                _make_assignment("SMN1", "AMBIGUOUS", 0.3),
            ],
        },
        variant_counts_by_locus={"SMN1": 5},
    )
    doc = json.loads(out.read_text())
    assert doc["loci"]["SMN1"]["status"] == "RESOLVED"


def test_summary_reports_empty_locus(tmp_path):
    out = tmp_path / "summary.json"
    write_summary(
        output_path=out,
        sample_name="HG002",
        reference="grch38",
        tech="ont",
        data_type="wgs",
        runtime_seconds=1.0,
        assignments_by_locus={"SMN1": []},
        variant_counts_by_locus={"SMN1": 0},
    )
    doc = json.loads(out.read_text())
    assert doc["loci"]["SMN1"]["status"] == "UNASSIGNED"
    assert doc["loci"]["SMN1"]["read_count"] == 0
