# LocusGuard

Caller-agnostic locus disambiguation engine for paralog and pseudogene regions — optimized for ONT long-read WGS, graceful on short-read.

## What it does

LocusGuard annotates variants in problematic genomic regions with the *true* locus of origin plus a confidence score, supported by multi-evidence scoring (PSV match, haplotype consistency, MAPQ pattern, soft-clip, coverage ratio) and reports homozygous-deletion status per paralog in the summary JSON. It does **not** call variants; it augments the output of your existing caller (Clair3, DeepVariant, etc.).

Each annotated variant carries:

- `TRUE_LOCUS` — assigned locus
- `LOCUS_CONF` — confidence in `[0, 1]`
- `LOCUS_STATUS` — `RESOLVED` | `PROBABLE` | `AMBIGUOUS` | `UNASSIGNED`
- `LOCUS_EVIDENCE` — per-source decomposition
- `LOCUS_KEY` — cross-output traceability key
- `GENE_CONVERSION_FLAG` — 1 if gene conversion suspected

## Status

**Phase 2.6-alt.** Supported: SMN1/SMN2 on GRCh38. ONT input with multi-evidence scoring (PSV match, haplotype consistency, MAPQ pattern, soft-clip, coverage ratio), haplotype clustering, gene-conversion detection, and **homozygous-deletion detection**. Short-read input is accepted but evidence adapters that require long reads are gracefully skipped.

Absolute copy-number quantitation was removed in this phase — benchmarking on HG002 R10.4.1 Dorado sup and olgu1 (R9/R10) showed the depth-based estimator is sensitive to aligner primary/secondary mapping decisions, yielding chemistry- and pipeline-dependent bias that cannot be closed with a static calibration factor. For absolute CN (SMA carrier screening etc.), use an orthogonal method (MLPA, ddPCR).

Not yet supported: WES mode (Phase 2.7), GBA/PMS2 (Phase 3), BAM output (Phase 3), packaging for Bioconda/Docker/nf-core (Phase 4).

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

- `patient.lg.vcf.gz` — annotated VCF with per-variant `LOCUS_*` INFO fields
- `patient.lg.summary.json` — per-locus status + `deletion_status` per locus
- `patient.lg.manifest.json` — versions, config hashes, command, runtime
- `patient.lg.report.html` — human-readable report with Deletion Status panel

## Library use

```python
from pathlib import Path
from locusguard.api import Annotator
from locusguard.config import load_config
from locusguard.deletion import classify_deletion

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

for locus_id, assignments in result.assignments_by_locus.items():
    print(f"{locus_id}: deletion_status = {classify_deletion(assignments)}")
```

## License

MIT
