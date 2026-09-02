import re
from pathlib import Path

FEATURE_MAPPER_PATH = (
    Path(__file__).parent.parent / "src" / "cfh" / "mapping" / "feature_mapper.py"
)

_LITERAL_PATTERN = re.compile(r'"P15056"|"NM_004333"|"BRAF"')


def test_feature_mapper_has_no_gene_specific_literals():
    matches = _LITERAL_PATTERN.findall(FEATURE_MAPPER_PATH.read_text())
    assert len(matches) == 0
