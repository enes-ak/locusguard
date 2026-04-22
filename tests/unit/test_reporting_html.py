from locusguard.reporting.html_report import write_html_report
from locusguard.types import Assignment, EvidenceScore, HaplotypeCluster


def _a(locus: str, status: str, conf: float) -> Assignment:
    return Assignment(
        read_id="r",
        assigned_locus=locus,
        confidence=conf,
        status=status,
        evidence_scores=[
            EvidenceScore(source="psv_match", normalized=conf, raw={}, available=True),
        ],
        locus_key=f"{locus}:abc",
        flags=set(),
    )


def test_html_contains_per_locus_cards(tmp_path):
    out = tmp_path / "report.html"
    write_html_report(
        output_path=out,
        sample_name="olgu1",
        reference="grch38",
        tech="ont",
        data_type="wgs",
        runtime_seconds=15.0,
        locusguard_version="0.2.0",
        assignments_by_locus={
            "SMN1": [_a("SMN1", "RESOLVED", 0.95)],
            "SMN2": [_a("SMN2", "RESOLVED", 0.90)],
        },
        clusters_by_locus={
            "SMN1": [HaplotypeCluster(
                hap_id="H1", supporting_reads=["r1", "r2"],
                psv_pattern={"P1": "C"}, assigned_locus="SMN1",
                confidence=0.95, notes=[],
            )],
            "SMN2": [],
        },
        variant_counts_by_locus={"SMN1": 1, "SMN2": 0},
        gene_conv_flags_by_locus={"SMN1": False, "SMN2": False},
        warnings=["Phase 1 scope"],
        degradations=[],
    )
    html = out.read_text()
    assert "<html" in html.lower()
    assert "olgu1" in html
    assert "SMN1" in html
    assert "SMN2" in html
    assert "RESOLVED" in html


def test_html_shows_gene_conversion_alert(tmp_path):
    out = tmp_path / "report.html"
    write_html_report(
        output_path=out,
        sample_name="s",
        reference="grch38",
        tech="ont",
        data_type="wgs",
        runtime_seconds=1.0,
        locusguard_version="0.2.0",
        assignments_by_locus={"SMN1": [_a("SMN1", "PROBABLE", 0.65)]},
        clusters_by_locus={"SMN1": [HaplotypeCluster(
            hap_id="H1", supporting_reads=["r1", "r2", "r3"],
            psv_pattern={"P1": "C", "P2": "G"}, assigned_locus="SMN1",
            confidence=0.5, notes=["gene_conversion_suspected", "hotspot_match:exon_7_conversion"],
        )]},
        variant_counts_by_locus={"SMN1": 1},
        gene_conv_flags_by_locus={"SMN1": True},
        warnings=[],
        degradations=[],
    )
    html = out.read_text()
    # The alert section should mention gene conversion and the hotspot
    assert "gene conversion" in html.lower() or "gene_conversion" in html.lower()
    assert "exon_7_conversion" in html


def test_html_renders_empty_locus(tmp_path):
    out = tmp_path / "empty.html"
    write_html_report(
        output_path=out,
        sample_name="s",
        reference="grch38",
        tech="ont",
        data_type="wgs",
        runtime_seconds=0.5,
        locusguard_version="0.2.0",
        assignments_by_locus={"SMN1": []},
        clusters_by_locus={"SMN1": []},
        variant_counts_by_locus={"SMN1": 0},
        gene_conv_flags_by_locus={"SMN1": False},
        warnings=[],
        degradations=[],
    )
    html = out.read_text()
    # Empty locus renders as UNASSIGNED status, no crash
    assert "SMN1" in html


def test_html_report_shows_deletion_status_present(tmp_path):
    """Present loci render a green PRESENT badge in the Deletion Status panel."""
    from locusguard.reporting.html_report import write_html_report
    from locusguard.types import Assignment

    assignments = [
        Assignment(
            read_id=f"r{i}", assigned_locus="SMN1", confidence=0.9,
            status="RESOLVED", evidence_scores=[], locus_key="SMN1:abc",
        )
        for i in range(12)
    ]
    out = tmp_path / "r.html"
    write_html_report(
        output_path=out, sample_name="x", reference="grch38",
        tech="ont", data_type="wgs", runtime_seconds=0.1,
        locusguard_version="test",
        assignments_by_locus={"SMN1": assignments},
        clusters_by_locus={"SMN1": []},
        variant_counts_by_locus={"SMN1": 0},
        gene_conv_flags_by_locus={"SMN1": False},
        warnings=[],
        degradations=[],
    )
    html = out.read_text()
    assert "Deletion Status" in html
    assert "PRESENT" in html


def test_html_report_shows_deletion_status_homozygous(tmp_path):
    from locusguard.reporting.html_report import write_html_report
    from locusguard.types import Assignment

    assignments = [
        Assignment(
            read_id=f"r{i}", assigned_locus=None, confidence=0.2,
            status="UNASSIGNED", evidence_scores=[], locus_key="SMN1:abc",
        )
        for i in range(15)
    ]
    out = tmp_path / "r.html"
    write_html_report(
        output_path=out, sample_name="x", reference="grch38",
        tech="ont", data_type="wgs", runtime_seconds=0.1,
        locusguard_version="test",
        assignments_by_locus={"SMN1": assignments},
        clusters_by_locus={"SMN1": []},
        variant_counts_by_locus={"SMN1": 0},
        gene_conv_flags_by_locus={"SMN1": False},
        warnings=[],
        degradations=[],
    )
    html = out.read_text()
    assert "HOMOZYGOUS_DELETION" in html
