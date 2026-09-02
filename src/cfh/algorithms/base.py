"""Plugin interface for hotspot-detection algorithms."""

from __future__ import annotations

from abc import ABC, abstractmethod

from cfh.genes.registry import GeneConfig
from cfh.model.algorithm_result import AlgorithmResult
from cfh.model.fusion_event import FusionEvent
from cfh.model.fusion_feature import FusionFeature


class Algorithm(ABC):
    """Base class every hotspot-detection algorithm plugin must implement."""

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
