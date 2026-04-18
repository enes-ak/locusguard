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
