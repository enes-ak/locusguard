"""Unit tests for capture bed parsing + PSV intersection."""
from __future__ import annotations

from pathlib import Path

import pytest

from locusguard.capture_bed import (
    CaptureBedError,
    CaptureRegion,
    PsvCoverage,
    compute_psv_coverage,
    load_capture_bed,
    position_in_capture,
)
from locusguard.config.schema import LocusConfig  # noqa: E402

# --- load_capture_bed ---

def test_load_bed3_basic(tmp_path: Path) -> None:
    p = tmp_path / "t.bed"
    p.write_text("chr5\t100\t200\nchr5\t300\t400\n")
    regions = load_capture_bed(p)
    assert regions == [CaptureRegion("chr5", 100, 200), CaptureRegion("chr5", 300, 400)]


def test_load_ignores_extra_columns(tmp_path: Path) -> None:
    p = tmp_path / "t.bed"
    p.write_text("chr5\t100\t200\tname1\t0\t+\nchr5\t300\t400\tname2\t0\t-\n")
    regions = load_capture_bed(p)
    assert regions == [CaptureRegion("chr5", 100, 200), CaptureRegion("chr5", 300, 400)]


def test_load_skips_track_browser_comment_blank(tmp_path: Path) -> None:
    p = tmp_path / "t.bed"
    p.write_text(
        "track name=example\n"
        "browser position chr5\n"
        "# a comment\n"
        "\n"
        "chr5\t100\t200\n"
    )
    regions = load_capture_bed(p)
    assert regions == [CaptureRegion("chr5", 100, 200)]


def test_load_malformed_row_raises(tmp_path: Path) -> None:
    p = tmp_path / "t.bed"
    p.write_text("chr5\tnotanumber\t200\n")
    with pytest.raises(CaptureBedError):
        load_capture_bed(p)


def test_load_missing_file_raises(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        load_capture_bed(tmp_path / "does_not_exist.bed")


def test_load_too_few_columns_raises(tmp_path: Path) -> None:
    p = tmp_path / "t.bed"
    p.write_text("chr5\t100\n")
    with pytest.raises(CaptureBedError):
        load_capture_bed(p)


# --- position_in_capture: boundary coverage ---

def test_position_inside() -> None:
    regions = [CaptureRegion("chr5", 100, 200)]
    assert position_in_capture("chr5", 150, regions) is True


def test_position_outside() -> None:
    regions = [CaptureRegion("chr5", 100, 200)]
    assert position_in_capture("chr5", 50, regions) is False
    assert position_in_capture("chr5", 250, regions) is False


def test_position_at_start_plus_one_is_covered() -> None:
    # 1-based pos == start+1 means 0-based pos == start, which IS in [start, end)
    regions = [CaptureRegion("chr5", 100, 200)]
    assert position_in_capture("chr5", 101, regions) is True


def test_position_at_end_is_covered() -> None:
    # 1-based pos == end means 0-based pos == end-1, which IS in [start, end)
    regions = [CaptureRegion("chr5", 100, 200)]
    assert position_in_capture("chr5", 200, regions) is True


def test_position_at_start_is_missing() -> None:
    # 1-based pos == start means 0-based pos == start-1, which is BEFORE [start, end)
    regions = [CaptureRegion("chr5", 100, 200)]
    assert position_in_capture("chr5", 100, regions) is False


def test_position_at_end_plus_one_is_missing() -> None:
    # 1-based pos == end+1 means 0-based pos == end, which is PAST [start, end)
    regions = [CaptureRegion("chr5", 100, 200)]
    assert position_in_capture("chr5", 201, regions) is False


def test_position_wrong_chrom() -> None:
    regions = [CaptureRegion("chr5", 100, 200)]
    assert position_in_capture("chr1", 150, regions) is False


def test_position_empty_regions() -> None:
    assert position_in_capture("chr5", 150, []) is False


# --- compute_psv_coverage ---

def _make_locus_config(psvs: list[tuple[str, int]]) -> LocusConfig:
    """Build a minimal LocusConfig with the given (name, pos_1based) PSVs."""
    return LocusConfig.model_validate({
        "schema_version": "1.0",
        "locusguard_compat": ">=0.1.0,<1.0.0",
        "locus": {
            "id": "TEST",
            "name": "TEST",
            "gene_family": "TEST",
            "paralogs": [],
        },
        "reference": "grch38",
        "coordinates": {
            "primary": {"chrom": "chr5", "start": 1, "end": 1000},
            "paralogs": {},
        },
        "psvs": [
            {"name": n, "chrom": "chr5", "pos": p, "alleles": {"TEST": "A"}}
            for n, p in psvs
        ],
        "evidence_weights": {
            "default": {
                "psv_match": 0.4,
                "haplotype_consistency": 0.25,
                "mapq_pattern": 0.1,
                "softclip": 0.05,
                "unique_kmer": 0.1,
                "coverage_ratio": 0.1,
            },
            "profile_overrides": {},
        },
        "confidence_thresholds": {"resolved": 0.8, "probable": 0.5},
    })


def test_compute_all_covered() -> None:
    cfg = _make_locus_config([("psv1", 150), ("psv2", 180)])
    regions = [CaptureRegion("chr5", 100, 200)]
    cov = compute_psv_coverage(cfg, regions)
    assert cov == PsvCoverage(covered=["psv1", "psv2"], missing=[], fraction_covered=1.0)


def test_compute_all_missing() -> None:
    cfg = _make_locus_config([("psv1", 50), ("psv2", 300)])
    regions = [CaptureRegion("chr5", 100, 200)]
    cov = compute_psv_coverage(cfg, regions)
    assert cov == PsvCoverage(covered=[], missing=["psv1", "psv2"], fraction_covered=0.0)


def test_compute_mixed() -> None:
    cfg = _make_locus_config([("inside", 150), ("outside", 300)])
    regions = [CaptureRegion("chr5", 100, 200)]
    cov = compute_psv_coverage(cfg, regions)
    assert cov.covered == ["inside"]
    assert cov.missing == ["outside"]
    assert cov.fraction_covered == 0.5


def test_compute_empty_psvs() -> None:
    cfg = _make_locus_config([])
    regions = [CaptureRegion("chr5", 100, 200)]
    cov = compute_psv_coverage(cfg, regions)
    assert cov == PsvCoverage(covered=[], missing=[], fraction_covered=0.0)


def test_compute_empty_regions() -> None:
    cfg = _make_locus_config([("psv1", 150)])
    cov = compute_psv_coverage(cfg, [])
    assert cov == PsvCoverage(covered=[], missing=["psv1"], fraction_covered=0.0)


def test_load_inverted_interval_raises(tmp_path: Path) -> None:
    """start > end is rejected as a malformed BED row."""
    p = tmp_path / "t.bed"
    p.write_text("chr5\t5000\t100\n")
    with pytest.raises(CaptureBedError, match="start .* > end"):
        load_capture_bed(p)
