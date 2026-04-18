import os
from pathlib import Path

import pytest
import yaml

from locusguard.io.reference import (
    ReferenceNotFoundError,
    resolve_reference_fasta,
)


def test_resolve_explicit_flag_wins(tmp_path, monkeypatch):
    fasta = tmp_path / "explicit.fa"
    fasta.write_text(">chr1\nACGT\n")
    env = tmp_path / "env.fa"
    env.write_text(">chr1\nTTTT\n")
    monkeypatch.setenv("LOCUSGUARD_GRCH38_FASTA", str(env))

    result = resolve_reference_fasta(
        reference="grch38",
        explicit_path=fasta,
        user_config_path=None,
    )
    assert result == fasta


def test_resolve_env_var_when_no_flag(tmp_path, monkeypatch):
    fasta = tmp_path / "env.fa"
    fasta.write_text(">chr1\nACGT\n")
    monkeypatch.setenv("LOCUSGUARD_GRCH38_FASTA", str(fasta))

    result = resolve_reference_fasta(
        reference="grch38",
        explicit_path=None,
        user_config_path=None,
    )
    assert result == fasta


def test_resolve_user_config_when_no_env_no_flag(tmp_path, monkeypatch):
    monkeypatch.delenv("LOCUSGUARD_GRCH38_FASTA", raising=False)
    fasta = tmp_path / "user.fa"
    fasta.write_text(">chr1\nACGT\n")
    user_cfg = tmp_path / "config.yaml"
    user_cfg.write_text(yaml.safe_dump({"references": {"grch38": str(fasta)}}))

    result = resolve_reference_fasta(
        reference="grch38",
        explicit_path=None,
        user_config_path=user_cfg,
    )
    assert result == fasta


def test_resolve_raises_when_none_found(tmp_path, monkeypatch):
    monkeypatch.delenv("LOCUSGUARD_GRCH38_FASTA", raising=False)
    with pytest.raises(ReferenceNotFoundError, match="grch38"):
        resolve_reference_fasta(
            reference="grch38",
            explicit_path=None,
            user_config_path=None,
        )


def test_resolve_raises_when_path_does_not_exist(tmp_path, monkeypatch):
    ghost = tmp_path / "does_not_exist.fa"
    monkeypatch.setenv("LOCUSGUARD_GRCH38_FASTA", str(ghost))
    with pytest.raises(ReferenceNotFoundError, match="does not exist"):
        resolve_reference_fasta(
            reference="grch38",
            explicit_path=None,
            user_config_path=None,
        )
