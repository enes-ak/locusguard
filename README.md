# LocusGuard

Caller-agnostic locus disambiguation engine for paralog and pseudogene regions — optimized for ONT long-read WGS, graceful on short-read.

## What it does

LocusGuard annotates variants in problematic genomic regions (paralog genes, pseudogenes, segmental duplications) with the *true* locus of origin plus a confidence score, supported by multi-evidence scoring: PSV match, haplotype consistency, MAPQ pattern, soft-clip pattern, and gene-conversion detection. It does **not** call variants; it augments the output of your existing caller (Clair3, DeepVariant, etc.).

Each annotated variant carries new INFO fields:

- `TRUE_LOCUS` — assigned locus (e.g. `SMN1`)
- `LOCUS_CONF` — confidence in `[0, 1]`
- `LOCUS_STATUS` — `RESOLVED` | `PROBABLE` | `AMBIGUOUS` | `UNASSIGNED`
- `LOCUS_EVIDENCE` — per-source decomposition
- `LOCUS_KEY` — cross-output traceability key
- `GENE_CONVERSION_FLAG` — 1 if gene conversion suspected at this locus

## Status

**Phase 2 (A-core).** Supported: SMN1/SMN2 on GRCh38. ONT + short-read inputs. Multi-evidence scoring with haplotype clustering and gene-conversion detection.

Not yet supported: GBA/PMS2 (Phase 3), BAM output (Phase 3), CN estimation, WES mode (Phase 2.5), packaging for Bioconda/Docker/nf-core (Phase 4).

## Install

```bash
pip install git+https://github.com/enes-ak/locusguard.git
```

Or from source (editable, for development):

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

- `patient.lg.vcf.gz` — annotated VCF
- `patient.lg.summary.json` — per-locus status, gene conversion flags, read counts
- `patient.lg.manifest.json` — versions, config hashes, command, runtime, degradations
- `patient.lg.report.html` — human-readable report (with `--emit-report`)
- `patient.lg.assignments.tsv` (with `--emit-assignments`) — per-read detail
- `patient.lg.haplotypes.tsv` (with `--emit-haplotypes`) — per-cluster detail

## Library use

```python
from pathlib import Path
from locusguard.api import Annotator
from locusguard.config import load_config

configs = [
    load_config("/path/to/SMN1.yaml"),
    load_config("/path/to/SMN2.yaml"),
]
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

# Phase 2 additions:
for locus_id, clusters in result.haplotype_clusters_by_locus.items():
    print(f"{locus_id}: {len(clusters)} haplotype cluster(s)")
```

## License

MIT
