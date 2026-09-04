import re
from pathlib import Path

MAPPING_DIR = Path(__file__).parent.parent / "src" / "cfh" / "mapping"
FEATURE_MAPPER_PATH = MAPPING_DIR / "feature_mapper.py"
GENOME_NEXUS_SOURCE_PATH = MAPPING_DIR / "genome_nexus_source.py"

_LITERAL_PATTERN = re.compile(r'"P15056"|"NM_004333"|"BRAF"')
_LITERAL_PATTERN_CASE_INSENSITIVE = re.compile(r"braf|P15056|NM_004333", re.IGNORECASE)


def test_feature_mapper_has_no_gene_specific_literals():
    matches = _LITERAL_PATTERN.findall(FEATURE_MAPPER_PATH.read_text())
    assert len(matches) == 0


def test_genome_nexus_source_has_no_gene_specific_literals():
    matches = _LITERAL_PATTERN_CASE_INSENSITIVE.findall(GENOME_NEXUS_SOURCE_PATH.read_text())
    assert len(matches) == 0
