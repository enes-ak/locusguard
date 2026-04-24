"""Shared pytest fixtures: build tiny BAM/VCF/FASTA files on the fly."""
from __future__ import annotations

from pathlib import Path

import pysam
import pytest

# Synthetic reference: a single chromosome `chr5` long enough to contain
# SMN1 and SMN2 approximate coordinates used in tests. We compress the 71 Mb
# full chromosome down to a 20 kb synthetic substitute that only keeps the
# two "locus regions" we care about, at relative offsets.
#
# Layout:
#   chr5: 0 .. 19999
#     SMN2 region:  2000 .. 8000
#     SMN1 region:  12000 .. 18000
#
# Tests requiring real GRCh38 coordinates use a separate fixture that
# builds a sparse-layout FASTA keyed to config coords.


MINI_CHROM = "chr5"
MINI_LENGTH = 20_000


def _random_dna(length: int, seed: int) -> str:
    import random
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(length))


@pytest.fixture
def mini_fasta(tmp_path: Path) -> Path:
    """Write a synthetic FASTA with a single small chr5."""
    path = tmp_path / "mini.fa"
    seq = _random_dna(MINI_LENGTH, seed=42)
    with path.open("w") as fh:
        fh.write(f">{MINI_CHROM}\n")
        for i in range(0, len(seq), 70):
            fh.write(seq[i : i + 70] + "\n")
    pysam.faidx(str(path))
    return path


@pytest.fixture
def smn_like_fasta(tmp_path: Path) -> Path:
    """A FASTA where two 3-kb windows mimic SMN1 and SMN2 regions.

    Both windows share >99% identity with one differing base at the
    'c.840C>T' discriminator position. The SMN1 window carries C; SMN2 T.

    Coordinates used match a compact test layout:
      chr5:  SMN2 at 2000..5000  (PSV at 4000: T)
             SMN1 at 12000..15000 (PSV at 14000: C)
    """
    path = tmp_path / "smn_like.fa"
    shared = _random_dna(3000, seed=7)
    seq = ["N"] * MINI_LENGTH
    seq[2000:5000] = list(shared)
    seq[12000:15000] = list(shared)
    seq[4000] = "T"            # SMN2 PSV allele
    seq[14000] = "C"           # SMN1 PSV allele
    with path.open("w") as fh:
        fh.write(f">{MINI_CHROM}\n")
        buf = "".join(seq)
        for i in range(0, len(buf), 70):
            fh.write(buf[i : i + 70] + "\n")
    pysam.faidx(str(path))
    return path


@pytest.fixture
def smn_like_bam(tmp_path: Path, smn_like_fasta: Path) -> Path:
    """A BAM with 10 reads over SMN1 PSV region + 10 reads over SMN2 PSV region.

    Reads are long (2 kb) ONT-style and all have MAPQ=60.
    SMN1-region reads carry base C at pos 14000.
    SMN2-region reads carry base T at pos 4000.
    """
    bam_path = tmp_path / "mini.bam"
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": MINI_CHROM, "LN": MINI_LENGTH}],
    }
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam:
        # SMN1-region reads
        for i in range(10):
            read = pysam.AlignedSegment()
            read.query_name = f"smn1_read_{i}"
            read.query_sequence = _read_sequence_for(smn_like_fasta, "chr5", 13000, 15000)
            read.flag = 0
            read.reference_id = 0
            read.reference_start = 13000
            read.mapping_quality = 60
            read.cigar = [(0, 2000)]  # 2000M
            read.query_qualities = pysam.qualitystring_to_array("I" * 2000)
            bam.write(read)

        # SMN2-region reads
        for i in range(10):
            read = pysam.AlignedSegment()
            read.query_name = f"smn2_read_{i}"
            read.query_sequence = _read_sequence_for(smn_like_fasta, "chr5", 3000, 5000)
            read.flag = 0
            read.reference_id = 0
            read.reference_start = 3000
            read.mapping_quality = 60
            read.cigar = [(0, 2000)]
            read.query_qualities = pysam.qualitystring_to_array("I" * 2000)
            bam.write(read)

    pysam.sort("-o", str(bam_path), str(bam_path))
    pysam.index(str(bam_path))
    return bam_path


def _read_sequence_for(fasta: Path, chrom: str, start: int, end: int) -> str:
    with pysam.FastaFile(str(fasta)) as fa:
        return fa.fetch(chrom, start, end).upper()


@pytest.fixture
def smn1_deleted_bam(tmp_path: Path, smn_like_fasta: Path) -> Path:
    """A BAM simulating a homozygous SMN1 deletion.

    15 reads span the SMN1 PSV region (ref pos ~14001) but carry the SMN2
    base (T) at the PSV position — simulating the clinical case where SMN1
    is deleted and SMN2 reads have been misaligned onto the SMN1 locus by
    the aligner. 10 normal SMN2 reads cover the SMN2 PSV region.

    With 15 reads all carrying the wrong PSV base, the SMN1 LocusAssigner
    should UNASSIGN every one of them (PSV mismatch → no confident
    SMN1-evidence), and ``classify_deletion`` should return
    ``HOMOZYGOUS_DELETION``.
    """
    bam_path = tmp_path / "smn1_deleted.bam"
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": MINI_CHROM, "LN": MINI_LENGTH}],
    }
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam:
        # 15 reads in SMN1 region (ref 13001-15000) with SMN2 PSV base (T) at 14001
        for i in range(15):
            read = pysam.AlignedSegment()
            read.query_name = f"smn1_deleted_read_{i}"
            seq = _read_sequence_for(smn_like_fasta, "chr5", 13000, 15000)
            # Offset of PSV pos 14001 (1-based) within sequence that starts at
            # ref pos 13001 (1-based) is 14001 - 13001 = 1000.
            # Replace that base with T (the SMN2 allele).
            seq = seq[:1000] + "T" + seq[1001:]
            read.query_sequence = seq
            read.flag = 0
            read.reference_id = 0
            read.reference_start = 13000
            read.mapping_quality = 60
            read.cigar = [(0, 2000)]  # 2000M
            read.query_qualities = pysam.qualitystring_to_array("I" * 2000)
            bam.write(read)

        # 10 reads in SMN2 region (ref 3001-5000) — normal SMN2 coverage
        for i in range(10):
            read = pysam.AlignedSegment()
            read.query_name = f"smn2_read_{i}"
            read.query_sequence = _read_sequence_for(smn_like_fasta, "chr5", 3000, 5000)
            read.flag = 0
            read.reference_id = 0
            read.reference_start = 3000
            read.mapping_quality = 60
            read.cigar = [(0, 2000)]
            read.query_qualities = pysam.qualitystring_to_array("I" * 2000)
            bam.write(read)

    pysam.sort("-o", str(bam_path), str(bam_path))
    pysam.index(str(bam_path))
    return bam_path


@pytest.fixture
def mini_vcf(tmp_path: Path) -> Path:
    """A gzipped, tabix-indexed VCF with two variants - one in SMN1 region, one in SMN2."""
    vcf_text = """##fileformat=VCFv4.2
##contig=<ID=chr5,length=20000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1
chr5\t3950\t.\tC\tT\t40\tPASS\t.\tGT\t0/1
chr5\t13950\t.\tA\tG\t45\tPASS\t.\tGT\t0/1
"""
    raw = tmp_path / "mini.vcf"
    raw.write_text(vcf_text)
    pysam.tabix_compress(str(raw), str(raw) + ".gz", force=True)
    pysam.tabix_index(str(raw) + ".gz", preset="vcf", force=True)
    return Path(str(raw) + ".gz")


@pytest.fixture
def multi_psv_fasta(tmp_path: Path) -> Path:
    """A FASTA where SMN1 and SMN2 each have 3 PSVs at known offsets.

    Layout (chr5):
      SMN1-like: pos 12000..15000
        PSV1 at 13000 = C (SMN1)
        PSV2 at 13500 = A (SMN1)
        PSV3 at 14000 = T (SMN1)
      SMN2-like: pos 2000..5000
        PSV1 at 3000 = T (SMN2)
        PSV2 at 3500 = G (SMN2)
        PSV3 at 4000 = A (SMN2)
    """
    path = tmp_path / "multi_psv.fa"
    seq = list(_random_dna(MINI_LENGTH, seed=101))
    # SMN1 PSVs
    seq[13000 - 1] = "C"
    seq[13500 - 1] = "A"
    seq[14000 - 1] = "T"
    # SMN2 PSVs
    seq[3000 - 1] = "T"
    seq[3500 - 1] = "G"
    seq[4000 - 1] = "A"
    with path.open("w") as fh:
        fh.write(f">{MINI_CHROM}\n")
        buf = "".join(seq)
        for i in range(0, len(buf), 70):
            fh.write(buf[i : i + 70] + "\n")
    pysam.faidx(str(path))
    return path


@pytest.fixture
def multi_psv_bam(tmp_path: Path, multi_psv_fasta: Path) -> Path:
    """A BAM with reads covering all 3 PSVs in each region."""
    bam_path = tmp_path / "multi_psv.bam"
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": MINI_CHROM, "LN": MINI_LENGTH}],
        "RG": [{"ID": "rg1", "PL": "ONT", "SM": "s1"}],
    }

    def _write_read(bam, name: str, start: int, length: int, fasta: Path):
        with pysam.FastaFile(str(fasta)) as fa:
            seq = fa.fetch(MINI_CHROM, start, start + length).upper()
        read = pysam.AlignedSegment()
        read.query_name = name
        read.query_sequence = seq
        read.flag = 0
        read.reference_id = 0
        read.reference_start = start
        read.mapping_quality = 60
        read.cigar = [(0, length)]
        read.query_qualities = pysam.qualitystring_to_array("I" * length)
        read.set_tag("RG", "rg1")
        bam.write(read)

    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam:
        for i in range(10):
            _write_read(bam, f"smn1_read_{i}", 12500, 2000, multi_psv_fasta)
        for i in range(10):
            _write_read(bam, f"smn2_read_{i}", 2500, 2000, multi_psv_fasta)

    pysam.sort("-o", str(bam_path), str(bam_path))
    pysam.index(str(bam_path))
    return bam_path


@pytest.fixture
def gene_conversion_bam(tmp_path: Path, multi_psv_fasta: Path) -> Path:
    """A BAM with 5 reads showing a mixed (gene-conversion-like) PSV pattern.

    Reads span the SMN1 region but carry:
      PSV1 at 13000 = 'T'  (SMN2-like)
      PSV2 at 13500 = 'A'  (SMN1-like)
      PSV3 at 14000 = 'T'  (SMN1-like)
    Mixed pattern should trigger gene conversion detection.
    """
    bam_path = tmp_path / "gene_conv.bam"
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": MINI_CHROM, "LN": MINI_LENGTH}],
        "RG": [{"ID": "rg1", "PL": "ONT", "SM": "s1"}],
    }
    with pysam.FastaFile(str(multi_psv_fasta)) as fa:
        base_seq = list(fa.fetch(MINI_CHROM, 12500, 14500).upper())
    # Flip PSV1 (1-based pos 13000 -> 0-based slice index 499) to T (SMN2-like)
    base_seq[499] = "T"
    synth_seq = "".join(base_seq)

    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam:
        for i in range(5):
            read = pysam.AlignedSegment()
            read.query_name = f"gc_read_{i}"
            read.query_sequence = synth_seq
            read.flag = 0
            read.reference_id = 0
            read.reference_start = 12500
            read.mapping_quality = 60
            read.cigar = [(0, 2000)]
            read.query_qualities = pysam.qualitystring_to_array("I" * 2000)
            read.set_tag("RG", "rg1")
            bam.write(read)
    pysam.sort("-o", str(bam_path), str(bam_path))
    pysam.index(str(bam_path))
    return bam_path
