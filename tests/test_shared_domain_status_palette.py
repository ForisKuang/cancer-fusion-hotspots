"""Guards against the domain-retention-status lollipop track
(``cfh.real_benchmark._domain_track_svg``) and the fusion-transcript
schematic (``cfh.reporting.fusion_schematic``) drifting apart on what a
color means.

Both renderers must draw their domain-retention-status colors from
:mod:`cfh.reporting.palette`. These tests fail (rather than silently
passing) if either renderer starts using a hardcoded hex literal again
instead of importing from that shared module.
"""

from __future__ import annotations

import inspect

from cfh import real_benchmark
from cfh.reporting import fusion_schematic, palette


def test_lollipop_track_imports_status_colors_from_shared_palette():
    source = inspect.getsource(real_benchmark)
    assert "from cfh.reporting.palette import" in source
    assert real_benchmark.RETAINED_COLOR is palette.RETAINED_COLOR
    assert real_benchmark.TRUNCATED_COLOR is palette.TRUNCATED_COLOR
    assert real_benchmark.LOST_COLOR is palette.LOST_COLOR
    assert real_benchmark.BREAKPOINT_COLOR is palette.BREAKPOINT_COLOR


def test_fusion_schematic_imports_status_colors_from_shared_palette():
    source = inspect.getsource(fusion_schematic)
    assert "from cfh.reporting.palette import" in source
    assert fusion_schematic.RETAINED_COLOR is palette.RETAINED_COLOR
    assert fusion_schematic.TRUNCATED_COLOR is palette.TRUNCATED_COLOR
    assert fusion_schematic.BREAKPOINT_COLOR is palette.BREAKPOINT_COLOR


def test_lollipop_status_dot_source_has_no_hardcoded_domain_status_hex_literals():
    """Neither renderer's own module source should hardcode the shared
    palette's hex values as string literals -- they must always come
    through the ``cfh.reporting.palette`` import, so a future edit that
    reintroduces a literal (even one that happens to match today) is
    caught."""
    for hex_value in (
        palette.RETAINED_COLOR,
        palette.TRUNCATED_COLOR,
        palette.LOST_COLOR,
        palette.BREAKPOINT_COLOR,
    ):
        assert f'"{hex_value}"' not in inspect.getsource(real_benchmark._domain_track_svg)


def test_fusion_schematic_source_has_no_hardcoded_domain_status_hex_literals():
    for hex_value in (
        palette.RETAINED_COLOR,
        palette.TRUNCATED_COLOR,
        palette.BREAKPOINT_COLOR,
    ):
        assert f'"{hex_value}"' not in inspect.getsource(fusion_schematic)
