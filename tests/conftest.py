"""Shared pytest fixtures and CLI options."""

from __future__ import annotations

import pytest


def pytest_addoption(parser: pytest.Parser) -> None:
    parser.addoption(
        "--depmap-path",
        action="store",
        default=None,
        help="Path to DepMap expression matrix (CSV/feather) for parity tests.",
    )


@pytest.fixture
def depmap_path(request: pytest.FixtureRequest) -> str | None:
    return request.config.getoption("--depmap-path")
