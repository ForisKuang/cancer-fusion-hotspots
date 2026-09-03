import pytest

from cfh.genes.registry import GeneConfig, load_gene_config


def test_loads_braf_yaml_into_validated_gene_config():
    config = load_gene_config("braf")
    assert isinstance(config, GeneConfig)
    assert config.canonical_transcript_id == "NM_004333"
    assert config.protein_id == "P15056"
    assert "kinase" in config.key_domains[0].name.lower()


def test_loads_ret_yaml_with_live_genome_nexus_identifiers():
    config = load_gene_config("ret")

    assert isinstance(config, GeneConfig)
    assert config.gene_symbol == "RET"
    assert config.canonical_transcript_id == "NM_020975"
    assert config.protein_id == "P07949"
    assert config.entrez_gene_id == 5979
    assert config.key_domains[0].accession == "PF07714"
    assert config.benchmark_reference.fusion_count == 7


def test_unknown_gene_raises():
    with pytest.raises(FileNotFoundError):
        load_gene_config("not_a_real_gene")
