from locusguard.io.vcf import VcfReader, VcfWriter


def test_vcf_reader_iterates_variants(mini_vcf):
    reader = VcfReader(mini_vcf)
    variants = list(reader.iter_variants())
    assert len(variants) == 2
    assert variants[0].CHROM == "chr5"
    assert variants[0].POS in (3950, 13950)  # cyvcf2 POS is 1-based


def test_vcf_reader_samples(mini_vcf):
    reader = VcfReader(mini_vcf)
    assert reader.samples == ["SAMPLE1"]


def test_vcf_reader_region_query(mini_vcf):
    reader = VcfReader(mini_vcf)
    in_smn1 = list(reader.fetch("chr5", 12000, 16000))
    assert len(in_smn1) == 1
    assert in_smn1[0].POS == 13950


def test_vcf_writer_round_trips_variants_with_info(mini_vcf, tmp_path):
    reader = VcfReader(mini_vcf)
    out_path = tmp_path / "out.vcf.gz"
    writer = VcfWriter(
        path=out_path,
        template_reader=reader,
        extra_info_fields=[
            ("TRUE_LOCUS", "1", "String", "Assigned locus"),
            ("LOCUS_CONF", "1", "Float", "Confidence [0,1]"),
        ],
    )
    for variant in reader.iter_variants():
        writer.write_annotated(variant, {"TRUE_LOCUS": "SMN1", "LOCUS_CONF": 0.87})
    writer.close()

    roundtrip = VcfReader(out_path)
    for variant in roundtrip.iter_variants():
        assert variant.INFO.get("TRUE_LOCUS") == "SMN1"
        assert abs(float(variant.INFO.get("LOCUS_CONF")) - 0.87) < 1e-6
