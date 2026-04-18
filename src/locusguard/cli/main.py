"""`locusguard` CLI entry point."""
from __future__ import annotations

import typer

from locusguard import __version__
from locusguard.cli.annotate import annotate

app = typer.Typer(
    name="locusguard",
    help="Caller-agnostic locus disambiguation engine for paralog and pseudogene regions.",
    no_args_is_help=True,
    add_completion=False,
)

app.command()(annotate)


@app.callback(invoke_without_command=True)
def _main(
    ctx: typer.Context,
    version: bool = typer.Option(
        False, "--version", "-V", help="Print version and exit.",
        is_eager=True,
    ),
) -> None:
    if version:
        typer.echo(__version__)
        raise typer.Exit()
    if ctx.invoked_subcommand is None:
        typer.echo(ctx.get_help())
        raise typer.Exit()


if __name__ == "__main__":
    app()
