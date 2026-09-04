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


@pytest.mark.parametrize(
    ("gene", "transcript", "protein", "entrez"),
    [
        ("alk", "NM_004304", "Q9UM73", 238),
        ("ntrk1", "NM_002529", "P04629", 4914),
    ],
)
def test_loads_alk_and_ntrk1_yaml_with_genome_nexus_canonical_identifiers(
    gene, transcript, protein, entrez
):
    """Curated IDs are the live Genome Nexus canonical-transcript values.

    Both targets use the catalytic protein-kinase Pfam family returned by
    that endpoint; its live residue bounds are deliberately resolved by the
    production Genome Nexus mapping path rather than copied into YAML.
    """
    config = load_gene_config(gene)

    assert config.canonical_transcript_id == transcript
    assert config.protein_id == protein
    assert config.entrez_gene_id == entrez
    assert config.key_domains[0].accession == "PF07714"
    assert config.key_domains[0].source == "genome_nexus"


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
