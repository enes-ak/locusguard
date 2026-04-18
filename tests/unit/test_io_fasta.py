import pytest

from locusguard.io.fasta import FastaReader


def test_fasta_reader_fetches_base_at_position(smn_like_fasta):
    reader = FastaReader(smn_like_fasta)
    # PSV at pos 14000 in SMN1 window = 'C'
    assert reader.base_at("chr5", 14000) == "C"
    # PSV at pos 4000 in SMN2 window = 'T'
    assert reader.base_at("chr5", 4000) == "T"


def test_fasta_reader_returns_region_sequence(smn_like_fasta):
    reader = FastaReader(smn_like_fasta)
    seq = reader.fetch("chr5", 2000, 2010)
    assert len(seq) == 10
    assert set(seq).issubset({"A", "C", "G", "T"})


def test_fasta_reader_list_chromosomes(smn_like_fasta):
    reader = FastaReader(smn_like_fasta)
    assert "chr5" in reader.chromosomes()


def test_fasta_reader_raises_on_unknown_chrom(smn_like_fasta):
    reader = FastaReader(smn_like_fasta)
    with pytest.raises(KeyError):
        reader.base_at("chr99", 0)
