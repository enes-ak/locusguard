# LocusGuard

Caller-agnostic locus disambiguation engine for paralog and pseudogene regions — optimized for ONT long-read WGS, graceful on short-read.

## What it does

LocusGuard annotates variants in problematic genomic regions (paralog genes, pseudogenes, segmental duplications) with the *true* locus of origin plus a confidence score. It does **not** call variants; it augments the output of your existing caller (Clair3, DeepVariant, etc.).

Each annotated variant carries new INFO fields:

- `TRUE_LOCUS` — assigned locus (e.g. `SMN1`)
- `LOCUS_CONF` — confidence in `[0, 1]`
- `LOCUS_STATUS` — `RESOLVED` | `PROBABLE` | `AMBIGUOUS` | `UNASSIGNED`
- `LOCUS_EVIDENCE` — per-source decomposition
- `LOCUS_KEY` — cross-output traceability key

## Status

**Phase 1 MVP.** Supported: SMN1/SMN2 on GRCh38. ONT + short-read inputs. PSV-based assignment.

Not yet supported: GBA/PMS2 (Phase 4), haplotype clustering (Phase 2), BAM output (Phase 3), HTML report (Phase 3), benchmark validation (Phase 5).

## Install

```bash
pip install locusguard  # once published to PyPI
# or, for development:
pip install -e ".[dev]"
```

## Quickstart

```bash
# Assumes GRCh38 FASTA available at $LOCUSGUARD_GRCH38_FASTA
export LOCUSGUARD_GRCH38_FASTA=/refs/grch38.primary.fa

locusguard annotate \
  --bam patient.bam \
  --vcf patient.vcf.gz \
  --reference grch38 \
  --tech ont \
  --data-type wgs \
  --locus SMN1,SMN2 \
  --output patient.lg.vcf.gz
```

Outputs:

- `patient.lg.vcf.gz` — annotated VCF (INFO fields added)
- `patient.lg.summary.json` — per-locus status and counts
- `patient.lg.manifest.json` — versions, config hashes, command, runtime

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
    summary_path=Path("patient.lg.summary.json"),
    manifest_path=Path("patient.lg.manifest.json"),
)
```

## Design

See [`docs/superpowers/specs/2026-04-18-locusguard-design.md`](docs/superpowers/specs/2026-04-18-locusguard-design.md) for the design rationale, architecture, and roadmap.

## License

MIT
