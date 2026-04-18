import pytest

from locusguard.assigner import LocusAssigner, build_analyzed_reads
from locusguard.config.schema import LocusConfig
from locusguard.io.bam import BamReader
from locusguard.io.fasta import FastaReader

# Reuse the fixture-compatible config from the PSV test, but with coords
# matching the smn_like_fasta / smn_like_bam fixtures.

_FIXTURE_SMN1_CFG = {
    "schema_version": "1.0",
    "locusguard_compat": ">=0.1.0,<1.0.0",
    "locus": {
        "id": "SMN1",
        "name": "SMN1",
        "gene_family": "SMN",
        "paralogs": ["SMN2"],
    },
    "reference": "grch38",
    "coordinates": {
        "primary": {"chrom": "chr5", "start": 12000, "end": 15000},
        "paralogs": {"SMN2": {"chrom": "chr5", "start": 2000, "end": 5000}},
    },
    "psvs": [
        {
            "name": "c.840C>T",
            "chrom": "chr5",
            "pos": 14000,
            "alleles": {"SMN1": "C", "SMN2": "T"},
        },
    ],
    "evidence_weights": {
        "default": {
            "psv_match": 0.40,
            "haplotype_consistency": 0.25,
            "mapq_pattern": 0.10,
            "softclip": 0.05,
            "unique_kmer": 0.10,
            "coverage_ratio": 0.10,
        },
        "profile_overrides": {},
    },
    "confidence_thresholds": {"resolved": 0.80, "probable": 0.50},
}

_FIXTURE_SMN2_CFG = {
    **_FIXTURE_SMN1_CFG,
    "locus": {
        "id": "SMN2",
        "name": "SMN2",
        "gene_family": "SMN",
        "paralogs": ["SMN1"],
    },
    "coordinates": {
        "primary": {"chrom": "chr5", "start": 2000, "end": 5000},
        "paralogs": {"SMN1": {"chrom": "chr5", "start": 12000, "end": 15000}},
    },
    "psvs": [
        {
            "name": "c.840C>T",
            "chrom": "chr5",
            "pos": 4000,
            "alleles": {"SMN1": "C", "SMN2": "T"},
        },
    ],
}


@pytest.fixture
def smn1_cfg() -> LocusConfig:
    return LocusConfig.model_validate(_FIXTURE_SMN1_CFG)


@pytest.fixture
def smn2_cfg() -> LocusConfig:
    return LocusConfig.model_validate(_FIXTURE_SMN2_CFG)


def test_build_analyzed_reads_extracts_psv_bases(smn_like_bam, smn_like_fasta, smn1_cfg):
    with BamReader(smn_like_bam) as bam, FastaReader(smn_like_fasta) as fa:
        reads = list(build_analyzed_reads(bam, fa, smn1_cfg, profile_name=None))
    assert len(reads) == 10
    for r in reads:
        psv = r.psv_observations["c.840C>T"]
        assert psv.reach is True
        # Reads fall in SMN1 region → observed base is C per fixture
        assert psv.base == "C"


def test_locus_assigner_assigns_smn1_reads_correctly(
    smn_like_bam, smn_like_fasta, smn1_cfg,
):
    assigner = LocusAssigner(smn1_cfg, profile_name=None)
    with BamReader(smn_like_bam) as bam, FastaReader(smn_like_fasta) as fa:
        assignments = assigner.assign(bam, fa)
    # 10 SMN1-region reads + 10 SMN2-region reads = 20 reads total pulled in
    # SMN1 reads have PSV=C → match → high conf; SMN2 reads have PSV=T → mismatch → low conf
    smn1_reads = [a for a in assignments if a.read_id.startswith("smn1_read_")]
    smn2_reads = [a for a in assignments if a.read_id.startswith("smn2_read_")]
    assert all(a.status == "RESOLVED" for a in smn1_reads)
    assert all(a.confidence >= 0.80 for a in smn1_reads)
    assert all(a.status in ("AMBIGUOUS", "UNASSIGNED") for a in smn2_reads)


def test_locus_key_is_stable(smn_like_bam, smn_like_fasta, smn1_cfg):
    assigner = LocusAssigner(smn1_cfg, profile_name=None)
    with BamReader(smn_like_bam) as bam, FastaReader(smn_like_fasta) as fa:
        first = assigner.assign(bam, fa)
    with BamReader(smn_like_bam) as bam, FastaReader(smn_like_fasta) as fa:
        second = assigner.assign(bam, fa)
    assert first[0].locus_key == second[0].locus_key
    assert first[0].locus_key.startswith("SMN1:")
