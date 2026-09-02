"""Built-in hotspot-detection algorithms."""

# Importing the module registers the built-in plugin with the shared registry.
from cfh.algorithms.domain_retention import DomainRetentionAlgorithm

__all__ = ["DomainRetentionAlgorithm"]
