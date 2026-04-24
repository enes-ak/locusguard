"""`locusguard annotate` — post-calling VCF annotation."""
from __future__ import annotations

from importlib.resources import files
from pathlib import Path
from typing import Annotated

import typer
from rich.console import Console

from locusguard.api import Annotator
from locusguard.config import load_config
from locusguard.config.schema import LocusConfig
from locusguard.io.reference import ReferenceNotFoundError, resolve_reference_fasta
from locusguard.preflight import PreflightError

console = Console()


def annotate(
    bam: Annotated[Path, typer.Option("--bam", help="Input BAM/CRAM with index.")],
    vcf: Annotated[Path, typer.Option("--vcf", help="Input VCF (bgzipped + tabix).")],
    output: Annotated[
        Path, typer.Option("--output", "-o", help="Output annotated VCF path (.vcf.gz).")
    ],
    reference: Annotated[
        str, typer.Option("--reference", help="Reference name (grch38).")
    ] = "grch38",
    reference_fasta: Annotated[
        Path | None,
        typer.Option("--reference-fasta", help="Path to reference FASTA (overrides env/config)."),
    ] = None,
    data_type: Annotated[
        str, typer.Option("--data-type", help="Data type (wgs | wes | panel).")
    ] = "wgs",
    capture_bed: Annotated[
        Path | None,
        typer.Option(
            "--capture-bed",
            help="BED file of capture regions (optional; valid with --data-type wes or panel).",
        ),
    ] = None,
    locus: Annotated[
        str | None,
        typer.Option(
            "--locus",
            help="Comma-separated locus IDs (SMN1,SMN2). Default: all bundled.",
        ),
    ] = None,
    config_paths: Annotated[
        list[Path] | None,
        typer.Option("--config", help="Additional user locus config YAML (repeatable)."),
    ] = None,
    emit_assignments: Annotated[
        bool, typer.Option("--emit-assignments", help="Write per-read assignments TSV.")
    ] = False,
    emit_haplotypes: Annotated[
        bool, typer.Option("--emit-haplotypes", help="Write per-cluster haplotypes TSV.")
    ] = False,
    emit_report: Annotated[
        bool, typer.Option("--emit-report", help="Write human-readable HTML report.")
    ] = False,
) -> None:
    """Post-calling VCF annotator (default LocusGuard mode).

    Adds per-variant INFO fields for locus disambiguation. Input VCF is read-only;
    output VCF is a new file with the original records plus TRUE_LOCUS, LOCUS_CONF,
    LOCUS_STATUS, LOCUS_EVIDENCE, LOCUS_KEY fields on variants inside configured
    locus regions.
    """
    valid_data_types = {"wgs", "wes", "panel"}
    if data_type not in valid_data_types:
        raise typer.BadParameter(
            f"--data-type must be one of: {', '.join(sorted(valid_data_types))}"
        )
    if data_type == "wgs" and capture_bed is not None:
        raise typer.BadParameter(
            "--capture-bed only valid with --data-type wes or panel"
        )
    if data_type in ("wes", "panel") and capture_bed is None:
        console.print(
            f"[yellow]notice:[/yellow] --data-type {data_type} without --capture-bed — "
            f"PSV coverage cannot be verified"
        )

    user_config_path = Path.home() / ".locusguard" / "config.yaml"
    try:
        fasta_path = resolve_reference_fasta(
            reference=reference,
            explicit_path=reference_fasta,
            user_config_path=user_config_path if user_config_path.exists() else None,
        )
    except ReferenceNotFoundError as e:
        console.print(f"[red]error:[/red] {e}")
        raise typer.Exit(1) from e

    configs = _load_configs(locus=locus, user_configs=config_paths or [])
    if not configs:
        console.print(
            "[red]error:[/red] no locus configs selected. "
            "Use --locus SMN1,SMN2 or --config path/to/yaml."
        )
        raise typer.Exit(1)

    annotator = Annotator(
        configs=configs,
        reference_fasta=fasta_path,
        data_type=data_type,
        capture_bed=capture_bed,
    )

    out_dir = output.parent
    out_dir.mkdir(parents=True, exist_ok=True)
    summary_path = output.with_suffix(".summary.json") if ".vcf" in output.name else None
    manifest_path = output.with_suffix(".manifest.json") if ".vcf" in output.name else None
    # Handle .vcf.gz two-suffix case
    if output.name.endswith(".vcf.gz"):
        base = output.name[: -len(".vcf.gz")]
        summary_path = out_dir / f"{base}.summary.json"
        manifest_path = out_dir / f"{base}.manifest.json"
    else:
        base = output.stem

    assignments_tsv_path = (
        out_dir / f"{base}.assignments.tsv" if emit_assignments else None
    )
    haplotypes_tsv_path = (
        out_dir / f"{base}.haplotypes.tsv" if emit_haplotypes else None
    )
    html_report_path = (
        out_dir / f"{base}.report.html" if emit_report else None
    )

    try:
        result = annotator.annotate_vcf(
            bam=bam,
            vcf_in=vcf,
            vcf_out=output,
            summary_path=summary_path,
            manifest_path=manifest_path,
            assignments_tsv_path=assignments_tsv_path,
            haplotypes_tsv_path=haplotypes_tsv_path,
            html_report_path=html_report_path,
        )
    except PreflightError as e:
        console.print(f"[red]preflight error:[/red] {e}")
        raise typer.Exit(1) from e

    console.print(
        f"[green]done:[/green] annotated {result.variants_annotated} of "
        f"{result.variants_total} variants across {len(configs)} loci. "
        f"Output: {output}"
    )


def _load_configs(locus: str | None, user_configs: list[Path]) -> list[LocusConfig]:
    configs: list[LocusConfig] = []

    # Bundled configs
    bundled_dir = files("locusguard").joinpath("configs/grch38/SMN")
    bundled_ids = {"SMN1", "SMN2"}
    wanted_ids = (
        set(locus.split(",")) if locus is not None else bundled_ids
    )
    for locus_id in sorted(wanted_ids & bundled_ids):
        path = bundled_dir.joinpath(f"{locus_id}.yaml")
        configs.append(load_config(str(path)))

    # User configs
    for p in user_configs:
        configs.append(load_config(p))

    return configs
