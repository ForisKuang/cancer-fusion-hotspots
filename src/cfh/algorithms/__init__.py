"""Algorithm plugin package.

Importing this package registers every built-in algorithm plugin (each
module below self-registers via the ``@register`` decorator on import), so
``registry.get(...)`` / ``registry.list_algorithms()`` work without callers
having to know which module defines a given algorithm name.
"""

from cfh.algorithms import confidence_stats  # noqa: F401
