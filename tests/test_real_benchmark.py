import json
from unittest.mock import MagicMock

from cfh.mapping.genome_nexus_source import GenomeNexusClient
from cfh.real_benchmark import analyze_structural_variant_calls, write_outputs


def test_real_benchmark_pipeline_writes_tsv_json_and_markdown(
    tmp_path,
    genome_nexus_canonical_transcript_fixture_path,
):
    client = MagicMock(spec=GenomeNexusClient)
    client.fetch_canonical_transcript.return_value = json.loads(
        genome_nexus_canonical_transcript_fixture_path.read_text()
    )
    calls = [
        {
            "sampleId": "SAMPLE-1",
            "site1HugoSymbol": "KIAA1549",
            "site2HugoSymbol": "BRAF",
            "site2Position": 140493152,
            "site2EffectOnFrame": "NA",
            "connectionType": "3to3",
            "eventInfo": "Protein Fusion: in frame  {KIAA1549:BRAF}",
        },
        {
            "sampleId": "SAMPLE-2",
            "site1HugoSymbol": "BRAF",
            "site2HugoSymbol": "AGK",
            "site1Position": 140493152,
            "site2EffectOnFrame": "NA",
            "connectionType": "5to5",
            "eventInfo": "Protein Fusion: out of frame  {BRAF:AGK}",
        },
    ]

    run = analyze_structural_variant_calls(
        calls,
        "BRAF",
        "msk_impact_50k_2026",
        genome_nexus_client=client,
        n_permutations=5,
    )
    paths = write_outputs(run, tmp_path, output_stem="benchmark")

    assert run.summary["total_fusions"] == 2
    assert run.summary["in_frame_count"] == 1
    assert paths["tsv"].read_text().splitlines()[0].startswith("event_id\tsample_id")
    payload = json.loads(paths["json"].read_text())
    assert payload["summary"]["domain_accession"] == "PF07714"
    assert len(payload["events"]) == 2
    report = paths["markdown"].read_text()
    assert "does **not** reproduce" in report
    assert "PF07714 (458-712 aa)" in report
