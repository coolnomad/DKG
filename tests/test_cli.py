"""Tests for dkg.cli."""

from __future__ import annotations

import json
from unittest.mock import patch

import pytest

from dkg.cli import main


def test_help_exits_zero():
    with pytest.raises(SystemExit) as exc_info:
        main(["--help"])
    assert exc_info.value.code == 0


def test_version_exits_zero(capsys):
    with pytest.raises(SystemExit) as exc_info:
        main(["--version"])
    assert exc_info.value.code == 0
    captured = capsys.readouterr()
    assert "dkg" in captured.out


def test_dispatch_xy(tmp_path):
    with patch("dkg.modes.xy.run") as mock_run:
        rc = main(
            [
                "--mode",
                "xy",
                "--x-matrix",
                "X.csv",
                "--y-matrix",
                "Y.csv",
                "--output-dir",
                str(tmp_path),
            ]
        )
    assert rc == 0
    mock_run.assert_called_once()
    config = mock_run.call_args[0][0]
    assert config.mode == "xy"
    assert config.x_matrix_path == "X.csv"
    assert config.y_matrix_path == "Y.csv"
    assert config.output_dir == str(tmp_path)


def test_dispatch_xx(tmp_path):
    with patch("dkg.modes.xx.run") as mock_run:
        rc = main(["--mode", "xx", "--x-matrix", "X.csv", "--output-dir", str(tmp_path)])
    assert rc == 0
    mock_run.assert_called_once()
    config = mock_run.call_args[0][0]
    assert config.mode == "xx"


def test_dispatch_pair(tmp_path):
    with patch("dkg.modes.pair.run") as mock_run:
        rc = main(
            [
                "--mode",
                "pair",
                "--x-matrix",
                "X.csv",
                "--y-matrix",
                "Y.csv",
                "--pair-x",
                "GENE_A",
                "--pair-y",
                "GENE_B",
                "--output-dir",
                str(tmp_path),
            ]
        )
    assert rc == 0
    mock_run.assert_called_once()
    config = mock_run.call_args[0][0]
    assert config.pair_x == "GENE_A"
    assert config.pair_y == "GENE_B"


def test_config_json_base_values(tmp_path):
    cfg_file = tmp_path / "cfg.json"
    cfg_file.write_text(json.dumps({"n_jobs": 4, "top_k": 500}))
    with patch("dkg.modes.xy.run") as mock_run:
        rc = main(
            [
                "--config-json",
                str(cfg_file),
                "--x-matrix",
                "X.csv",
                "--y-matrix",
                "Y.csv",
                "--output-dir",
                str(tmp_path),
            ]
        )
    assert rc == 0
    config = mock_run.call_args[0][0]
    assert config.n_jobs == 4
    assert config.top_k == 500


def test_config_json_cli_overrides(tmp_path):
    cfg_file = tmp_path / "cfg.json"
    cfg_file.write_text(json.dumps({"top_k": 500}))
    with patch("dkg.modes.xy.run") as mock_run:
        rc = main(
            [
                "--config-json",
                str(cfg_file),
                "--top-k",
                "999",
                "--x-matrix",
                "X.csv",
                "--y-matrix",
                "Y.csv",
                "--output-dir",
                str(tmp_path),
            ]
        )
    assert rc == 0
    config = mock_run.call_args[0][0]
    assert config.top_k == 999


def test_tier1_threshold_sets_both(tmp_path):
    with patch("dkg.modes.xy.run") as mock_run:
        rc = main(
            [
                "--tier1-threshold",
                "0.35",
                "--x-matrix",
                "X.csv",
                "--y-matrix",
                "Y.csv",
                "--output-dir",
                str(tmp_path),
            ]
        )
    assert rc == 0
    config = mock_run.call_args[0][0]
    assert config.tier1_pearson_threshold == pytest.approx(0.35)
    assert config.tier1_spearman_threshold == pytest.approx(0.35)


def test_runner_error_returns_1(tmp_path, capsys):
    with patch("dkg.modes.xy.run", side_effect=RuntimeError("boom")):
        rc = main(["--x-matrix", "X.csv", "--y-matrix", "Y.csv", "--output-dir", str(tmp_path)])
    assert rc == 1
    assert "boom" in capsys.readouterr().err


def test_bad_config_json_returns_1(tmp_path, capsys):
    bad_file = tmp_path / "bad.json"
    bad_file.write_text("not json {{{")
    rc = main(["--config-json", str(bad_file), "--output-dir", str(tmp_path)])
    assert rc == 1
    assert "error" in capsys.readouterr().err.lower()


def test_missing_config_json_returns_1(tmp_path, capsys):
    rc = main(["--config-json", str(tmp_path / "nonexistent.json"), "--output-dir", str(tmp_path)])
    assert rc == 1


def test_invalid_mode_exits_one(capsys):
    with pytest.raises(SystemExit) as exc_info:
        main(["--mode", "invalid"])
    assert exc_info.value.code == 1
    assert "error" in capsys.readouterr().err.lower()


def test_default_mode_is_xy(tmp_path):
    with patch("dkg.modes.xy.run") as mock_run:
        rc = main(["--x-matrix", "X.csv", "--y-matrix", "Y.csv", "--output-dir", str(tmp_path)])
    assert rc == 0
    config = mock_run.call_args[0][0]
    assert config.mode == "xy"
