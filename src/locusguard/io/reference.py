"""Resolve reference FASTA path via flag, env var, or user config."""
from __future__ import annotations

import os
from pathlib import Path

import yaml

_ENV_VAR_BY_REFERENCE = {
    "grch38": "LOCUSGUARD_GRCH38_FASTA",
}


class ReferenceNotFoundError(FileNotFoundError):
    """Raised when a reference FASTA cannot be resolved."""


def resolve_reference_fasta(
    reference: str,
    explicit_path: Path | None,
    user_config_path: Path | None,
) -> Path:
    """Resolve a reference FASTA path using the documented priority chain.

    Priority: 1) explicit_path, 2) env var, 3) user config file.

    Raises ReferenceNotFoundError if none yields an existing path.
    """
    candidate = _pick_candidate(reference, explicit_path, user_config_path)
    if candidate is None:
        env_var = _ENV_VAR_BY_REFERENCE.get(reference, f"LOCUSGUARD_{reference.upper()}_FASTA")
        raise ReferenceNotFoundError(
            f"Reference FASTA for '{reference}' not provided. "
            f"Set one of: --reference-fasta /path/to.fa, "
            f"{env_var} env var, or add a 'references.{reference}: /path/to.fa' "
            f"entry to ~/.locusguard/config.yaml."
        )
    if not candidate.exists():
        raise ReferenceNotFoundError(
            f"Reference FASTA path does not exist: {candidate}"
        )
    return candidate


def _pick_candidate(
    reference: str,
    explicit_path: Path | None,
    user_config_path: Path | None,
) -> Path | None:
    if explicit_path is not None:
        return Path(explicit_path)

    env_var = _ENV_VAR_BY_REFERENCE.get(reference)
    if env_var:
        env_value = os.environ.get(env_var)
        if env_value:
            return Path(env_value)

    if user_config_path is not None and user_config_path.exists():
        with user_config_path.open("r", encoding="utf-8") as fh:
            doc = yaml.safe_load(fh) or {}
        refs = doc.get("references", {})
        if reference in refs:
            return Path(refs[reference])

    return None
