from locusguard.depth import DepthStats, compute_region_depth
from locusguard.io.bam import BamReader


def test_depth_stats_fields():
    d = DepthStats(
        region_name="test",
        chrom="chr5",
        start=1000,
        end=2000,
        mean_depth=25.5,
        median_depth=26.0,
        length_bp=1000,
        reads_counted=30,
    )
    assert d.length_bp == 1000
    assert d.mean_depth == 25.5


def test_compute_region_depth_on_covered_region(smn_like_bam):
    with BamReader(smn_like_bam) as bam:
        # smn_like_bam has 10 reads (2000 bp each) starting at pos 13000
        stats = compute_region_depth(bam, "chr5", 13500, 14500, "smn1_core")
    assert stats.region_name == "smn1_core"
    assert stats.length_bp == 1001  # inclusive end
    assert stats.mean_depth > 5.0  # ~10x with some partial coverage
    assert stats.reads_counted >= 10


def test_compute_region_depth_on_empty_region(smn_like_bam):
    with BamReader(smn_like_bam) as bam:
        stats = compute_region_depth(bam, "chr5", 19000, 19500, "empty")
    assert stats.mean_depth == 0.0
    assert stats.reads_counted == 0


def test_compute_region_depth_on_unknown_chrom(smn_like_bam):
    import pytest
    with BamReader(smn_like_bam) as bam, pytest.raises((ValueError, KeyError)):
        compute_region_depth(bam, "chrUnknown", 1, 100, "test")
