from locusguard.cn import estimate_cn
from locusguard.config.schema import LocusConfig
from locusguard.depth.region import DepthStats


def _minimal_cfg(control_count: int = 2) -> LocusConfig:
    cfg = {
        "schema_version": "1.0",
        "locusguard_compat": ">=0.1.0,<1.0.0",
        "locus": {"id": "SMN1", "name": "SMN1", "gene_family": "SMN", "paralogs": ["SMN2"]},
        "reference": "grch38",
        "coordinates": {
            "primary": {"chrom": "chr5", "start": 100, "end": 500},
            "paralogs": {"SMN2": {"chrom": "chr5", "start": 1000, "end": 1400}},
        },
        "psvs": [
            {
                "name": "P1", "chrom": "chr5", "pos": 150,
                "alleles": {"SMN1": "C", "SMN2": "T"},
            }
        ],
        "evidence_weights": {
            "default": {
                "psv_match": 0.40, "haplotype_consistency": 0.25, "mapq_pattern": 0.10,
                "softclip": 0.05, "unique_kmer": 0.10, "coverage_ratio": 0.10,
            },
            "profile_overrides": {},
        },
        "confidence_thresholds": {"resolved": 0.80, "probable": 0.50},
        "control_regions": [
            {
                "name": f"ctrl_{i}", "chrom": "chr5",
                "start": 10000 + i * 1000, "end": 10999 + i * 1000,
            }
            for i in range(control_count)
        ],
    }
    return LocusConfig.model_validate(cfg)


def _depth(name: str, mean: float, length: int = 400) -> DepthStats:
    return DepthStats(
        region_name=name, chrom="chr5", start=1, end=length,
        mean_depth=mean, median_depth=mean, length_bp=length,
        reads_counted=int(mean * length / 2000),
    )


def test_cn_estimator_normal_diploid():
    cfg = _minimal_cfg()
    depths = {
        "SMN1": _depth("SMN1", 30.0),
        "SMN2": _depth("SMN2", 30.0),
        "ctrl_0": _depth("ctrl_0", 30.0),
        "ctrl_1": _depth("ctrl_1", 30.0),
    }
    est = estimate_cn(cfg, depths, tech="ont")
    assert est.status == "ok"
    # 2 * (30 / 30) = 2
    assert abs(est.absolute_cn - 2.0) < 0.01
    assert est.absolute_cn_rounded == 2
    assert abs(est.paralog_ratio - 1.0) < 0.01


def test_cn_estimator_deletion():
    cfg = _minimal_cfg()
    depths = {
        "SMN1": _depth("SMN1", 0.0),       # SMN1 deleted
        "SMN2": _depth("SMN2", 30.0),
        "ctrl_0": _depth("ctrl_0", 30.0),
        "ctrl_1": _depth("ctrl_1", 30.0),
    }
    est = estimate_cn(cfg, depths, tech="ont")
    assert est.status == "ok"
    assert est.absolute_cn_rounded == 0


def test_cn_estimator_insufficient_depth():
    cfg = _minimal_cfg()
    depths = {
        "SMN1": _depth("SMN1", 3.0),       # below 10x threshold
        "SMN2": _depth("SMN2", 30.0),
        "ctrl_0": _depth("ctrl_0", 30.0),
    }
    est = estimate_cn(cfg, depths, tech="ont")
    assert est.status == "insufficient_depth"
    assert est.absolute_cn is None


def test_cn_estimator_no_control_regions():
    cfg = _minimal_cfg(control_count=0)
    depths = {
        "SMN1": _depth("SMN1", 30.0),
        "SMN2": _depth("SMN2", 30.0),
    }
    est = estimate_cn(cfg, depths, tech="ont")
    assert est.status == "no_control_regions"
    assert est.method == "no_control_regions"


def test_cn_estimator_unsupported_tech():
    cfg = _minimal_cfg()
    depths = {
        "SMN1": _depth("SMN1", 30.0),
        "SMN2": _depth("SMN2", 30.0),
        "ctrl_0": _depth("ctrl_0", 30.0),
    }
    est = estimate_cn(cfg, depths, tech="short-read")
    assert est.status == "not_supported_for_tech"
    assert est.method == "not_supported_for_tech"


def test_cn_estimator_rounds_fractional():
    cfg = _minimal_cfg()
    depths = {
        "SMN1": _depth("SMN1", 45.0),       # 2 * 45/30 = 3
        "SMN2": _depth("SMN2", 15.0),       # 2 * 15/30 = 1
        "ctrl_0": _depth("ctrl_0", 30.0),
        "ctrl_1": _depth("ctrl_1", 30.0),
    }
    est = estimate_cn(cfg, depths, tech="ont")
    assert est.absolute_cn_rounded == 3
    assert abs(est.paralog_ratio - 3.0) < 0.01
