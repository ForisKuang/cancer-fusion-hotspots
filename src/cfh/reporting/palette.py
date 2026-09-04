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

import colorsys
import zlib

RETAINED_COLOR = "#2878b5"
"""A domain (or domain segment) that is fully retained."""

TRUNCATED_COLOR = "#f2a93b"
"""A domain (or domain segment) that is partially retained/truncated."""

LOST_COLOR = "#777777"
"""A domain that is fully lost."""

BREAKPOINT_COLOR = "#d62728"
"""Breakpoint marker; also used as the lollipop track's
reference-discrepancy outline color."""


def deterministic_color(label: str, *, lightness: float = 0.55, saturation: float = 0.55) -> str:
    """Deterministic, arbitrary-but-stable hex color for an arbitrary
    string label: the same label always gets the same color within a run
    and across runs (a pure function of the label), so a reader can
    visually track one entity across a diagram. Uses a CRC32 hash into
    hue space rather than Python's salted ``hash()``, which is randomized
    per-process and would make the same label render a different color on
    every regeneration.

    Shared by :func:`cfh.reporting.fusion_schematic.partner_color`
    (partner-gene coloring) and
    :func:`cfh.real_benchmark._domain_track_svg` (per-domain highlight
    coloring when a gene configures more than one key domain), so the two
    renderers can't drift apart on how a stable color is derived.
    """
    hue = (zlib.crc32(label.encode("utf-8")) % 360) / 360.0
    r, g, b = colorsys.hls_to_rgb(hue, lightness, saturation)
    return f"#{int(r * 255):02x}{int(g * 255):02x}{int(b * 255):02x}"
