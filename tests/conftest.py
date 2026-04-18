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
