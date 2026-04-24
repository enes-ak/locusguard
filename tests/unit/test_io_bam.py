from locusguard.io.bam import BamReader


def test_bam_reader_fetches_reads_in_region(smn_like_bam):
    reader = BamReader(smn_like_bam)
    reads = list(reader.fetch("chr5", 12000, 16000))
    assert len(reads) == 10
    assert all(r.query_name.startswith("smn1_read_") for r in reads)


def test_bam_reader_excludes_out_of_region_reads(smn_like_bam):
    reader = BamReader(smn_like_bam)
    reads = list(reader.fetch("chr5", 2000, 6000))
    assert len(reads) == 10
    assert all(r.query_name.startswith("smn2_read_") for r in reads)


def test_bam_reader_header_accessible(smn_like_bam):
    reader = BamReader(smn_like_bam)
    assert reader.is_sorted_by_coordinate is True
    assert "chr5" in reader.chromosomes()


def test_bam_reader_estimated_is_long_read_returns_true_for_long_reads(smn_like_bam):
    reader = BamReader(smn_like_bam)
    # Our fixture reads are 2000 bp → classified as long
    assert reader.estimated_is_long_read() is True


