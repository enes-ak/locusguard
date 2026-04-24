import pytest

from locusguard.preflight import PreflightError, run_preflight


def test_preflight_passes_on_valid_inputs(smn_like_bam, mini_vcf, smn_like_fasta):
    # Should not raise.
    run_preflight(bam=smn_like_bam, vcf=mini_vcf, fasta=smn_like_fasta)


def test_preflight_fails_on_missing_bam_index(smn_like_bam, mini_vcf, smn_like_fasta, tmp_path):
    # Copy BAM without index
    unindexed = tmp_path / "unindexed.bam"
    unindexed.write_bytes(smn_like_bam.read_bytes())
    with pytest.raises(PreflightError, match="index"):
        run_preflight(bam=unindexed, vcf=mini_vcf, fasta=smn_like_fasta)


def test_preflight_fails_on_chrom_mismatch(smn_like_bam, mini_vcf, tmp_path):
    # FASTA with different chromosome name
    bad_fasta = tmp_path / "bad.fa"
    bad_fasta.write_text(">chrOther\nACGTACGT\n")
    import pysam
    pysam.faidx(str(bad_fasta))
    with pytest.raises(PreflightError, match="chromosome"):
        run_preflight(bam=smn_like_bam, vcf=mini_vcf, fasta=bad_fasta)


def test_preflight_returns_none_without_capture_bed(smn_like_bam, smn_like_fasta, mini_vcf):
    from locusguard.preflight import run_preflight
    result = run_preflight(bam=smn_like_bam, vcf=mini_vcf, fasta=smn_like_fasta)
    assert result is None


def test_preflight_loads_valid_capture_bed(tmp_path, smn_like_bam, smn_like_fasta, mini_vcf):
    from locusguard.capture_bed import CaptureRegion
    from locusguard.preflight import run_preflight
    bed = tmp_path / "cap.bed"
    bed.write_text("chr5\t100\t200\nchr5\t300\t400\n")
    result = run_preflight(
        bam=smn_like_bam, vcf=mini_vcf, fasta=smn_like_fasta, capture_bed=bed,
    )
    assert result == [CaptureRegion("chr5", 100, 200), CaptureRegion("chr5", 300, 400)]


def test_preflight_missing_capture_bed_raises(tmp_path, smn_like_bam, smn_like_fasta, mini_vcf):
    import pytest

    from locusguard.preflight import PreflightError, run_preflight
    bed = tmp_path / "does_not_exist.bed"
    with pytest.raises(PreflightError, match="Capture bed file not found"):
        run_preflight(
            bam=smn_like_bam, vcf=mini_vcf, fasta=smn_like_fasta, capture_bed=bed,
        )


def test_preflight_malformed_capture_bed_raises(tmp_path, smn_like_bam, smn_like_fasta, mini_vcf):
    import pytest

    from locusguard.preflight import PreflightError, run_preflight
    bed = tmp_path / "bad.bed"
    bed.write_text("chr5\tnotanumber\t200\n")
    with pytest.raises(PreflightError, match="Failed to parse capture bed"):
        run_preflight(
            bam=smn_like_bam, vcf=mini_vcf, fasta=smn_like_fasta, capture_bed=bed,
        )
