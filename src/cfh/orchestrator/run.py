"""Run registered algorithms concurrently against a common input snapshot."""

from __future__ import annotations

import hashlib
import json
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from time import perf_counter
from typing import Any, Sequence

from cfh.algorithms.registry import get
from cfh.genes.registry import GeneConfig
from cfh.model.algorithm_result import AlgorithmResult
from cfh.model.fusion_event import FusionEvent
from cfh.model.fusion_feature import FusionFeature


def _input_fingerprint(events: Sequence[FusionEvent], features: Sequence[FusionFeature]) -> str:
    """Hash a canonical representation of the input tables."""
    payload = {
        "events": [event.model_dump(mode="json") for event in events],
        "features": [feature.model_dump(mode="json") for feature in features],
    }
    serialized = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(serialized.encode()).hexdigest()


def _params_for(
    name: str,
    params: dict[str, Any] | None,
    algorithm_names: Sequence[str] | None = None,
) -> dict[str, Any]:
    """Resolve optional per-algorithm parameters without sharing mutable state."""
    if not params:
        return {}
    named_params = params.get(name)
    if isinstance(named_params, dict):
        return dict(named_params)
    if algorithm_names and any(
        isinstance(params.get(algorithm), dict) for algorithm in algorithm_names
    ):
        return {}
    if all(isinstance(value, dict) for value in params.values()):
        return {}
    return dict(params)


def _dependencies_for(name: str) -> tuple[str, ...]:
    """Return the ``DEPENDS_ON`` names declared for a registered algorithm.

    An unregistered ``name`` has no dependencies here -- ``_run_one`` still
    surfaces the lookup failure as that algorithm's own failed result.
    """
    try:
        algorithm_cls = get(name)
    except KeyError:
        return ()
    return tuple(getattr(algorithm_cls, "DEPENDS_ON", ()) or ())


def _schedule_waves(algorithm_names: Sequence[str]) -> list[list[str]]:
    """Group distinct algorithm names into dependency-respecting waves via
    Kahn's algorithm (topological sort).

    An algorithm with no declared ``DEPENDS_ON`` (the default for every
    existing plugin) always lands in the first wave it is eligible for,
    preserving today's fully-concurrent dispatch for the common case. An
    algorithm that declares dependencies is deferred to the first wave in
    which every also-requested dependency has already completed, so their
    results can be injected into its params before it runs. A dependency
    that was not itself requested never blocks scheduling -- it is simply
    unavailable to the dependent algorithm.

    Raises ``ValueError`` naming the algorithms involved if a circular
    ``DEPENDS_ON`` dependency is found among the requested algorithms
    (e.g. A depends on B and B depends on A) -- a configuration error that
    can never make scheduling progress, so it is reported explicitly
    rather than silently degraded into one wave with unresolved,
    empty-injected dependencies.
    """
    remaining = list(dict.fromkeys(algorithm_names))  # de-duplicated, order-preserving
    scheduled: set[str] = set()
    waves: list[list[str]] = []
    while remaining:
        wave = [
            name
            for name in remaining
            if all(dep in scheduled or dep not in remaining for dep in _dependencies_for(name))
        ]
        if not wave:
            raise ValueError(
                "circular DEPENDS_ON dependency detected among algorithms: "
                f"{', '.join(sorted(remaining))}"
            )
        waves.append(wave)
        scheduled.update(wave)
        remaining = [name for name in remaining if name not in wave]
    return waves


def _json_safe_params(params: dict[str, Any]) -> dict[str, Any]:
    """Coerce a params dict to something ``AlgorithmResult.model_dump(mode="json")``
    can always serialize.

    A failed algorithm's synthetic result (below) echoes back its raw input
    ``params`` verbatim as ``Parameters`` for debuggability -- but a caller
    may have passed a live object in there (e.g. a shared
    ``GenomeNexusClient``, as ``cutpoint_detection``/``domain_retention``
    already accept) purely for the algorithm's own internal use, never
    intended to be serialized. Any value that is not JSON-safe is replaced
    with its ``repr()`` so serialization never breaks on a legitimately
    failed algorithm result.
    """

    def _safe(value: Any) -> Any:
        if value is None or isinstance(value, (str, int, float, bool)):
            return value
        if isinstance(value, dict):
            return {str(key): _safe(item) for key, item in value.items()}
        if isinstance(value, (list, tuple)):
            return [_safe(item) for item in value]
        return repr(value)

    return _safe(params)


def _run_one(
    name: str,
    events: list[FusionEvent],
    features: list[FusionFeature],
    gene_config: GeneConfig,
    params: dict[str, Any],
    input_fingerprint: str,
) -> AlgorithmResult:
    """Execute one plugin and convert exceptions to structured failures."""
    started_at = perf_counter()
    try:
        result = get(name)().run(events, features, gene_config, params)
        summary = dict(result.Summary or {})
        summary["Runtime_seconds"] = perf_counter() - started_at
        return result.model_copy(
            update={
                "Parameters": result.Parameters if result.Parameters is not None else params,
                "Summary": summary,
                "Warnings": result.Warnings or [],
                "Created_at": result.Created_at or datetime.now(timezone.utc),
                "Input_fingerprint": input_fingerprint,
            }
        )
    except Exception as exc:
        return AlgorithmResult(
            Algorithm=name,
            Parameters=_json_safe_params(params),
            Summary={"Runtime_seconds": perf_counter() - started_at},
            Tables={},
            Warnings=[f"Algorithm failed: {type(exc).__name__}: {exc}"],
            Created_at=datetime.now(timezone.utc),
            Input_fingerprint=input_fingerprint,
        )


def run_algorithms(
    algorithm_names: Sequence[str],
    events: list[FusionEvent],
    features: list[FusionFeature],
    gene_config: GeneConfig,
    params: dict[str, Any] | None = None,
    *,
    max_workers: int | None = None,
    extra_results: Sequence[AlgorithmResult] | None = None,
) -> list[AlgorithmResult]:
    """Run registered algorithms and return results in request order.

    ``params`` may contain per-algorithm dictionaries keyed by algorithm name.
    ``results_to_json`` serializes this return value as a combined JSON array.
    Exceptions from an individual plugin become that plugin's warning result,
    allowing independent algorithms in the same run to complete successfully.

    Algorithms are dispatched in dependency-respecting waves (see
    :func:`_schedule_waves`): the common case of no declared
    :attr:`~cfh.algorithms.base.Algorithm.DEPENDS_ON` runs everything in one
    fully concurrent wave, exactly as before. An algorithm that *does*
    declare dependencies runs in a later wave, once its also-requested
    dependencies have completed, with their results automatically merged
    into its ``params["algorithm_results"]`` (unless the caller already set
    that key explicitly). ``extra_results`` seeds this dependency pool with
    results computed outside this call (e.g. by a caller that pre-computed
    one algorithm separately before requesting the rest), without re-running
    or re-returning them.
    """
    input_fingerprint = _input_fingerprint(events, features)
    if not algorithm_names:
        return []

    indexed_names = list(enumerate(algorithm_names))
    waves = _schedule_waves(algorithm_names)

    results: list[AlgorithmResult | None] = [None] * len(algorithm_names)
    results_by_name: dict[str, AlgorithmResult] = {
        result.Algorithm: result for result in (extra_results or [])
    }

    for wave in waves:
        wave_entries = [(index, name) for index, name in indexed_names if name in wave]
        wave_params: dict[int, dict[str, Any]] = {}
        for index, name in wave_entries:
            algorithm_params = _params_for(name, params, algorithm_names)
            depends_on = _dependencies_for(name)
            if depends_on and "algorithm_results" not in algorithm_params:
                available = [results_by_name[dep] for dep in depends_on if dep in results_by_name]
                if available:
                    algorithm_params = {**algorithm_params, "algorithm_results": available}
            wave_params[index] = algorithm_params

        worker_count = max_workers or len(wave_entries)
        with ThreadPoolExecutor(max_workers=worker_count) as executor:
            futures = {
                executor.submit(
                    _run_one,
                    name,
                    [event.model_copy(deep=True) for event in events],
                    [feature.model_copy(deep=True) for feature in features],
                    gene_config.model_copy(deep=True),
                    wave_params[index],
                    input_fingerprint,
                ): index
                for index, name in wave_entries
            }
            for future in as_completed(futures):
                index = futures[future]
                result = future.result()
                results[index] = result
                results_by_name[result.Algorithm] = result

    return [result for result in results if result is not None]


def results_to_json(results: Sequence[AlgorithmResult]) -> str:
    """Serialize combined algorithm results into a JSON array for emission."""
    return json.dumps([result.model_dump(mode="json") for result in results])
