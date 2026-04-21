# LocusGuard

Caller-agnostic locus disambiguation engine for paralog and pseudogene regions — optimized for ONT long-read WGS, graceful on short-read.

## What it does

LocusGuard annotates variants in problematic genomic regions with the *true* locus of origin plus a confidence score, supported by multi-evidence scoring (PSV match, haplotype consistency, MAPQ, soft-clip, coverage ratio) and reports absolute copy number per paralog. It does **not** call variants; it augments the output of your existing caller (Clair3, DeepVariant, etc.).

Each annotated variant carries:

- `TRUE_LOCUS` — assigned locus
- `LOCUS_CONF` — confidence in `[0, 1]`
- `LOCUS_STATUS` — `RESOLVED` | `PROBABLE` | `AMBIGUOUS` | `UNASSIGNED`
- `LOCUS_EVIDENCE` — per-source decomposition
- `LOCUS_KEY` — cross-output traceability key
- `GENE_CONVERSION_FLAG` — 1 if gene conversion suspected
- `CN_CONTEXT` — per-locus absolute copy number

## Status

**Phase 2.5.** Supported: SMN1/SMN2 on GRCh38. ONT input for CN estimation; short-read graceful with CN disabled. Multi-evidence scoring with haplotype clustering and gene-conversion detection. Control-region-normalized copy number.

Not yet supported: WES mode (Phase 2.6), GBA/PMS2 (Phase 3), BAM output (Phase 3), packaging for Bioconda/Docker/nf-core (Phase 4).

## Install

```bash
pip install git+https://github.com/enes-ak/locusguard.git
```

Or from source:

```bash
git clone https://github.com/enes-ak/locusguard
cd locusguard
pip install -e ".[dev]"
```

## Quickstart

```bash
export LOCUSGUARD_GRCH38_FASTA=/refs/grch38.primary.fa

locusguard annotate \
  --bam patient.bam \
  --vcf patient.vcf.gz \
  --reference grch38 \
  --tech ont \
  --data-type wgs \
  --locus SMN1,SMN2 \
  --emit-report \
  --output patient.lg.vcf.gz
```

Outputs:

- `patient.lg.vcf.gz` — annotated VCF (INFO fields including CN_CONTEXT)
- `patient.lg.summary.json` — per-locus status + `cn_estimate` block
- `patient.lg.manifest.json` — versions, config hashes, command, runtime, CN method + controls used
- `patient.lg.report.html` — human-readable report with Copy Number panel

## Library use

```python
from pathlib import Path
from locusguard.api import Annotator
from locusguard.config import load_config

configs = [load_config(p) for p in ["SMN1.yaml", "SMN2.yaml"]]
annotator = Annotator(
    configs=configs,
    reference_fasta=Path("/refs/grch38.primary.fa"),
    tech="ont",
    data_type="wgs",
)
result = annotator.annotate_vcf(
    bam=Path("patient.bam"),
    vcf_in=Path("patient.vcf.gz"),
    vcf_out=Path("patient.lg.vcf.gz"),
    html_report_path=Path("patient.lg.report.html"),
)

for locus_id, cn in result.cn_by_locus.items():
    if cn.status == "ok":
        print(f"{locus_id}: CN = {cn.absolute_cn_rounded} (est: {cn.absolute_cn:.2f})")
```

## License

MIT
