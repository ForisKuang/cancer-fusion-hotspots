"""Registration mechanism for hotspot-detection algorithm plugins."""

from __future__ import annotations

from cfh.algorithms.base import Algorithm

_REGISTRY: dict[str, type[Algorithm]] = {}


def register(name: str):
    """Class decorator that registers an :class:`Algorithm` under ``name``."""

    def _decorator(cls: type[Algorithm]) -> type[Algorithm]:
        if not issubclass(cls, Algorithm):
            raise TypeError(f"{cls!r} must subclass Algorithm to be registered")
        _REGISTRY[name] = cls
        return cls

    return _decorator


def get(name: str) -> type[Algorithm]:
    """Look up a registered algorithm class by name."""
    try:
        return _REGISTRY[name]
    except KeyError as exc:
        raise KeyError(f"No algorithm registered under {name!r}") from exc


def list_algorithms() -> list[str]:
    """List all registered algorithm names."""
    return sorted(_REGISTRY)


def unregister(name: str) -> None:
    """Remove a registered algorithm (mainly useful for test isolation)."""
    _REGISTRY.pop(name, None)
