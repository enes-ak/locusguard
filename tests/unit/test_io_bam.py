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


def test_bam_reader_detects_short_reads_heuristically(smn_like_bam):
    reader = BamReader(smn_like_bam)
    # Our fixture reads are 2000 bp → classified as long
    assert reader.estimated_is_long_read() is True


def test_detect_tech_from_rg_pl_illumina(tmp_path):
    import pysam

    bam_path = tmp_path / "illumina.bam"
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr5", "LN": 20000}],
        "RG": [{"ID": "rg1", "PL": "ILLUMINA", "SM": "s1"}],
    }
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam:
        read = pysam.AlignedSegment()
        read.query_name = "r1"
        read.query_sequence = "A" * 150
        read.flag = 0
        read.reference_id = 0
        read.reference_start = 1000
        read.mapping_quality = 40
        read.cigar = [(0, 150)]
        read.query_qualities = pysam.qualitystring_to_array("I" * 150)
        read.set_tag("RG", "rg1")
        bam.write(read)
    pysam.sort("-o", str(bam_path), str(bam_path))
    pysam.index(str(bam_path))

    from locusguard.io.bam import BamReader
    with BamReader(bam_path) as reader:
        assert reader.detect_tech() == "short-read"


def test_detect_tech_from_rg_pl_ont(tmp_path):
    import pysam

    bam_path = tmp_path / "ont.bam"
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr5", "LN": 20000}],
        "RG": [{"ID": "rg1", "PL": "ONT", "SM": "s1"}],
    }
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam:
        read = pysam.AlignedSegment()
        read.query_name = "r1"
        read.query_sequence = "A" * 5000
        read.flag = 0
        read.reference_id = 0
        read.reference_start = 1000
        read.mapping_quality = 60
        read.cigar = [(0, 5000)]
        read.query_qualities = pysam.qualitystring_to_array("I" * 5000)
        read.set_tag("RG", "rg1")
        bam.write(read)
    pysam.sort("-o", str(bam_path), str(bam_path))
    pysam.index(str(bam_path))

    from locusguard.io.bam import BamReader
    with BamReader(bam_path) as reader:
        assert reader.detect_tech() == "ont"


def test_detect_tech_no_rg_fallback_to_length(smn_like_bam):
    from locusguard.io.bam import BamReader
    with BamReader(smn_like_bam) as reader:
        # Our fixture reads are 2000 bp → long-read
        assert reader.detect_tech() == "ont"
