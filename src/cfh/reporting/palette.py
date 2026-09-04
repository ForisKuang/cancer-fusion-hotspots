"""Single source of truth for the domain-retention-status SVG color
convention shared by every renderer that visualizes it.

Two renderers currently draw from this palette:
:func:`cfh.real_benchmark._domain_track_svg` (the per-event domain-retention
lollipop track) and :mod:`cfh.reporting.fusion_schematic` (the
fusion-transcript schematic). Both must import these constants rather than
hardcoding their own hex literals, so the two visualizations can never
silently drift apart on what a color means.
"""

from __future__ import annotations

RETAINED_COLOR = "#2878b5"
"""A domain (or domain segment) that is fully retained."""

TRUNCATED_COLOR = "#f2a93b"
"""A domain (or domain segment) that is partially retained/truncated."""

LOST_COLOR = "#777777"
"""A domain that is fully lost."""

BREAKPOINT_COLOR = "#d62728"
"""Breakpoint marker; also used as the lollipop track's
reference-discrepancy outline color."""
