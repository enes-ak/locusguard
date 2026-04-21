from locusguard.haplotype import cluster_reads
from locusguard.types import AnalyzedRead, PSVObs


def _read(read_id: str, **psvs_with_reach) -> AnalyzedRead:
    obs = {}
    for name, (base, reach) in psvs_with_reach.items():
        obs[name] = PSVObs(base=base, qual=30, reach=reach)
    return AnalyzedRead(
        read_id=read_id,
        aligned_chrom="chr5",
        aligned_pos=1000,
        psv_observations=obs,
        mapq=60,
        softclip_5p=0,
        softclip_3p=0,
        is_long_read=True,
        is_supplementary=False,
        original_mapq_zero=False,
    )


def test_identical_patterns_group_into_one_cluster():
    reads = [
        _read("r1", P1=("C", True), P2=("C", True), P3=("T", True)),
        _read("r2", P1=("C", True), P2=("C", True), P3=("T", True)),
    ]
    clusters = cluster_reads(reads, psv_names=["P1", "P2", "P3"])
    assert len(clusters) == 1
    assert clusters[0].psv_pattern == {"P1": "C", "P2": "C", "P3": "T"}
    assert set(clusters[0].supporting_reads) == {"r1", "r2"}


def test_different_patterns_form_separate_clusters():
    reads = [
        _read("r1", P1=("C", True), P2=("C", True)),
        _read("r2", P1=("T", True), P2=("T", True)),
        _read("r3", P1=("C", True), P2=("C", True)),
    ]
    clusters = cluster_reads(reads, psv_names=["P1", "P2"])
    assert len(clusters) == 2
    sizes = sorted(len(c.supporting_reads) for c in clusters)
    assert sizes == [1, 2]


def test_missing_reach_encoded_as_question_mark():
    reads = [
        _read("r1", P1=("C", True), P2=("?", False)),
    ]
    clusters = cluster_reads(reads, psv_names=["P1", "P2"])
    assert clusters[0].psv_pattern == {"P1": "C", "P2": "?"}


def test_cluster_ids_are_stable_and_sorted():
    reads = [
        _read("r1", P1=("T", True)),
        _read("r2", P1=("A", True)),
        _read("r3", P1=("C", True)),
    ]
    clusters = cluster_reads(reads, psv_names=["P1"])
    assert [c.hap_id for c in clusters] == ["H1", "H2", "H3"]
    patterns = [c.psv_pattern["P1"] for c in clusters]
    assert patterns == sorted(patterns)


def test_empty_reads_returns_empty():
    assert cluster_reads([], psv_names=["P1"]) == []
