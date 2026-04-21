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


def test_summary_includes_gene_conv_flag(tmp_path):
    import json

    from locusguard.reporting.summary import write_summary

    out = tmp_path / "summary.json"
    write_summary(
        output_path=out,
        sample_name="s",
        reference="grch38",
        tech="ont",
        data_type="wgs",
        runtime_seconds=1.0,
        assignments_by_locus={
            "SMN1": [_make_assignment("SMN1", "PROBABLE", 0.65)],
        },
        variant_counts_by_locus={"SMN1": 1},
        gene_conv_flags_by_locus={"SMN1": True},
        gene_conv_evidence_by_locus={"SMN1": "mixed PSV pattern at hotspot:exon_7_conversion"},
    )
    doc = json.loads(out.read_text())
    assert doc["loci"]["SMN1"]["gene_conv_flag"] is True
    assert "exon_7_conversion" in doc["loci"]["SMN1"]["gene_conv_evidence"]


def test_summary_includes_cn_estimate_block(tmp_path):
    import json

    from locusguard.reporting.summary import write_summary
    from locusguard.types import CnEstimate

    out = tmp_path / "summary.json"
    write_summary(
        output_path=out,
        sample_name="s",
        reference="grch38",
        tech="ont",
        data_type="wgs",
        runtime_seconds=1.0,
        assignments_by_locus={"SMN1": [_make_assignment("SMN1", "RESOLVED", 0.9)]},
        variant_counts_by_locus={"SMN1": 1},
        cn_by_locus={
            "SMN1": CnEstimate(
                locus_id="SMN1",
                absolute_cn=2.05,
                absolute_cn_rounded=2,
                paralog_ratio=2.0,
                cn_total_family=3.05,
                method="control_region_normalized",
                confidence=0.87,
                status="ok",
                notes=["control:chr5_A"],
            ),
        },
    )
    doc = json.loads(out.read_text())
    cn = doc["loci"]["SMN1"]["cn_estimate"]
    assert cn["absolute_cn"] == 2.05
    assert cn["absolute_cn_rounded"] == 2
    assert cn["method"] == "control_region_normalized"
    assert cn["status"] == "ok"


def test_summary_cn_estimate_skipped_when_no_cn(tmp_path):
    import json

    from locusguard.reporting.summary import write_summary

    out = tmp_path / "summary.json"
    write_summary(
        output_path=out,
        sample_name="s",
        reference="grch38",
        tech="ont",
        data_type="wgs",
        runtime_seconds=1.0,
        assignments_by_locus={"SMN1": [_make_assignment("SMN1", "RESOLVED", 0.9)]},
        variant_counts_by_locus={"SMN1": 1},
    )
    doc = json.loads(out.read_text())
    # CN block optional — absent when cn_by_locus not provided
    assert "cn_estimate" not in doc["loci"]["SMN1"]
