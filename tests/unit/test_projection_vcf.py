
from locusguard.io.vcf import VcfReader
from locusguard.projection.vcf import VcfProjector
from locusguard.types import Assignment, EvidenceScore


def _synthetic_assignments(locus: str) -> list[Assignment]:
    return [
        Assignment(
            read_id=f"r{i}",
            assigned_locus=locus,
            confidence=0.87,
            status="RESOLVED",
            evidence_scores=[
                EvidenceScore(
                    source="psv_match",
                    normalized=0.87,
                    raw={"matches": 2, "total_observations": 2},
                    available=True,
                ),
            ],
            locus_key=f"{locus}:a3f9c2",
            flags=set(),
        )
        for i in range(5)
    ]


def test_projector_writes_all_expected_info_fields(mini_vcf, tmp_path):
    out = tmp_path / "out.vcf.gz"
    projector = VcfProjector(
        input_vcf=mini_vcf,
        output_vcf=out,
        locus_regions=[
            ("SMN1", "chr5", 12000, 15000),
            ("SMN2", "chr5", 2000, 5000),
        ],
    )
    assignments_by_locus = {
        "SMN1": _synthetic_assignments("SMN1"),
        "SMN2": _synthetic_assignments("SMN2"),
    }
    projector.run(assignments_by_locus)

    reader = VcfReader(out)
    variants = list(reader.iter_variants())
    assert len(variants) == 2

    annotated_fields = ("TRUE_LOCUS", "LOCUS_CONF", "LOCUS_STATUS",
                        "LOCUS_EVIDENCE", "LOCUS_KEY")
    for v in variants:
        for f in annotated_fields:
            assert v.INFO.get(f) is not None, f"{f} missing on {v.CHROM}:{v.POS}"
        assert v.INFO.get("LOCUS_STATUS") == "RESOLVED"


def test_variant_outside_any_region_is_passed_through(mini_vcf, tmp_path):
    """Variants outside configured locus regions carry no LocusGuard INFO."""
    out = tmp_path / "out.vcf.gz"
    projector = VcfProjector(
        input_vcf=mini_vcf,
        output_vcf=out,
        locus_regions=[
            # Region that doesn't overlap either fixture variant
            ("SMN1", "chr5", 19000, 19500),
        ],
    )
    projector.run({"SMN1": _synthetic_assignments("SMN1")})

    reader = VcfReader(out)
    variants = list(reader.iter_variants())
    for v in variants:
        # INFO.get returns None for missing fields with cyvcf2
        assert v.INFO.get("TRUE_LOCUS") is None


def test_projection_emits_gene_conversion_flag(mini_vcf, tmp_path):
    out = tmp_path / "out.vcf.gz"
    assignments = [
        Assignment(
            read_id=f"r{i}",
            assigned_locus="SMN1",
            confidence=0.70,
            status="PROBABLE",
            evidence_scores=[],
            locus_key="SMN1:a3f9c2",
            flags={"gene_conversion_suspected"},
        )
        for i in range(5)
    ]
    projector = VcfProjector(
        input_vcf=mini_vcf,
        output_vcf=out,
        locus_regions=[("SMN1", "chr5", 12000, 15000)],
    )
    projector.run({"SMN1": assignments})

    reader = VcfReader(out)
    variants = list(reader.iter_variants())
    smn1_variant = next(v for v in variants if v.POS == 13950)
    assert int(smn1_variant.INFO.get("GENE_CONVERSION_FLAG")) == 1



