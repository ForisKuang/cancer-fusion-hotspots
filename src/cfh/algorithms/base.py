"""Plugin interface for hotspot-detection algorithms."""

from __future__ import annotations

from abc import ABC, abstractmethod

from cfh.genes.registry import GeneConfig
from cfh.model.algorithm_result import AlgorithmResult
from cfh.model.fusion_event import FusionEvent
from cfh.model.fusion_feature import FusionFeature


class Algorithm(ABC):
    """Base class every hotspot-detection algorithm plugin must implement."""

    DEPENDS_ON: tuple[str, ...] = ()
    """Names of other registered algorithms whose results this algorithm
    consumes. Empty by default -- the vast majority of algorithms are
    independent and the orchestrator dispatches them all concurrently, as
    before. An algorithm that declares a non-empty ``DEPENDS_ON`` is instead
    scheduled by :func:`cfh.orchestrator.run.run_algorithms` in a later wave,
    once every also-requested dependency has completed, with those
    dependencies' ``AlgorithmResult`` objects automatically injected into its
    own ``params["algorithm_results"]`` (unless the caller already supplied
    that key explicitly). A dependency that was not itself requested in the
    same call is simply unavailable, not a scheduling error -- it is up to
    the dependent algorithm to handle a missing/partial set gracefully.
    """

    @abstractmethod
    def run(
        self,
        events: list[FusionEvent],
        features: list[FusionFeature],
        gene_config: GeneConfig,
        params: dict,
    ) -> AlgorithmResult:
        """Run the algorithm and return a structured :class:`AlgorithmResult`."""
        raise NotImplementedError
