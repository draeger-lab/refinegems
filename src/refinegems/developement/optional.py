"""Helpers for optional runtime dependencies."""

from __future__ import annotations

import importlib
from types import ModuleType


class OptionalDependencyError(ImportError):
    """Raised when an optional feature dependency is not installed."""


def require_optional_dependency(
    import_name: str,
    *,
    package_name: str | None = None,
    extra: str | None = None,
    purpose: str | None = None,
) -> ModuleType:
    """Import an optional dependency or raise an actionable error."""
    try:
        return importlib.import_module(import_name)
    except ImportError as exc:
        install_name = package_name or import_name
        message = f"Optional dependency '{install_name}' is required"
        if purpose:
            message += f" to {purpose}"
        message += "."
        if extra:
            message += f" Install refineGEMs with the '{extra}' extra."
        raise OptionalDependencyError(message) from exc


def require_bioregistry(purpose: str | None = None) -> ModuleType:
    """Import bioregistry for annotation-related features."""
    return require_optional_dependency(
        "bioregistry",
        extra="bioregistry",
        purpose=purpose or "use annotation registry features",
    )
