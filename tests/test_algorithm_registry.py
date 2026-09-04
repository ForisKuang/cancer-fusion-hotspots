from cfh.algorithms import registry
from cfh.algorithms.base import Algorithm
from cfh.model.algorithm_result import AlgorithmResult


def test_register_and_retrieve_dummy_algorithm():
    @registry.register("dummy")
    class DummyAlgorithm(Algorithm):
        def run(self, events, features, gene_config, params) -> AlgorithmResult:
            return AlgorithmResult(Algorithm="dummy")

    try:
        assert registry.get("dummy") is DummyAlgorithm
        assert "dummy" in registry.list_algorithms()

        result = DummyAlgorithm().run([], [], None, {})
        assert result.Algorithm == "dummy"
    finally:
        registry.unregister("dummy")
