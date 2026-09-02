"""Built-in hotspot-detection algorithm plugins.

Importing this package registers every built-in algorithm plugin via each
module's ``@register`` decorator.
"""

from cfh.algorithms.confidence_stats import ConfidenceStatsAlgorithm
from cfh.algorithms.exon_retention import ExonRetentionAnalysis
from cfh.algorithms.frequency import FrequencyAnalysis

__all__ = [
    "ConfidenceStatsAlgorithm",
    "ExonRetentionAnalysis",
    "FrequencyAnalysis",
]
