import pytest

from cfh.genes.registry import GeneConfig, load_gene_config


def test_loads_braf_yaml_into_validated_gene_config():
    config = load_gene_config("braf")
    assert isinstance(config, GeneConfig)
    assert config.canonical_transcript_id == "NM_004333"
    assert config.protein_id == "P15056"
    assert "kinase" in config.key_domains[0].name.lower()


def test_unknown_gene_raises():
    with pytest.raises(FileNotFoundError):
        load_gene_config("not_a_real_gene")


def test_braf_disruption_required_domains_are_configured_with_real_pfam_accessions():
    """BRAF opts into the domain-disruption test with its real N-terminal
    autoinhibitory module (RAS-binding + cysteine-rich domains), sourced
    from the same live Genome Nexus data used for the kinase domain."""
    config = load_gene_config("braf")
    accessions = {domain.accession for domain in config.disruption_required_domains}
    assert accessions == {"PF02196", "PF00130"}
    assert all(domain.source == "genome_nexus" for domain in config.disruption_required_domains)


def test_disruption_required_domains_defaults_to_empty_list():
    """Opt-in field: a config that never mentions it must gracefully default
    to no-op, not error, the same as the existing key_domains pattern."""
    config = GeneConfig(gene_symbol="FAKE", canonical_transcript_id="NM_1", protein_id="P1")
    assert config.disruption_required_domains == []


def test_eml4_alk_and_a_synthetic_tmprss2_erg_pair_leave_disruption_domains_unset():
    """Negative controls: EML4-ALK (joint-partner mechanism) and TMPRSS2-ERG
    (promoter-swap/expression-driven mechanism) must not configure this
    field, so domain_disruption gracefully no-ops for both."""
    eml4_alk = load_gene_config("eml4-alk")
    assert eml4_alk.disruption_required_domains == []

    tmprss2_erg = GeneConfig(gene_pair=("TMPRSS2", "ERG"))
    assert tmprss2_erg.disruption_required_domains == []
