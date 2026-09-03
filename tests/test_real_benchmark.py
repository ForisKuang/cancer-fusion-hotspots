import json
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock

import pytest
import requests
from click.testing import CliRunner

from cfh import cli
from cfh.ingestion import cbioportal_api
from cfh.mapping.genome_nexus_source import GenomeNexusClient
from cfh.real_benchmark import (
    RealBenchmarkNetworkError,
    analyze_structural_variant_calls,
    run_real_benchmark,
    write_outputs,
)


def _genome_nexus_client(fixture_path):
    client = MagicMock(spec=GenomeNexusClient)
    client.fetch_canonical_transcript.return_value = json.loads(fixture_path.read_text())
    return client


def _call(
    sample_id,
    *,
    event_info="Protein Fusion: in frame  {KIAA1549:BRAF}",
    breakpoint=140493152,
    connection_type="3to3",
):
    return {
        "sampleId": sample_id,
        "site1HugoSymbol": "KIAA1549",
        "site2HugoSymbol": "BRAF",
        "site2Position": breakpoint,
        "site2EffectOnFrame": "NA",
        "connectionType": connection_type,
        "eventInfo": event_info,
    }


def test_real_benchmark_pipeline_writes_tsv_json_and_markdown(
    tmp_path,
    genome_nexus_canonical_transcript_fixture_path,
):
    client = _genome_nexus_client(genome_nexus_canonical_transcript_fixture_path)
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
    assert run.summary["mapped_fusions"] == 2
    assert run.summary["in_frame_count"] == 1
    assert run.warnings == []
    assert paths["tsv"].read_text().splitlines()[0].startswith("event_id\tsample_id")
    payload = json.loads(paths["json"].read_text())
    assert payload["summary"]["domain_accession"] == "PF07714"
    assert len(payload["events"]) == 2
    report = paths["markdown"].read_text()
    assert "does **not** reproduce" in report
    assert "PF07714 (458-712 aa)" in report


def test_malformed_fusion_rows_are_warned_and_skipped_without_losing_valid_rows(
    genome_nexus_canonical_transcript_fixture_path,
):
    client = _genome_nexus_client(genome_nexus_canonical_transcript_fixture_path)
    calls = [
        _call("VALID"),
        _call(
            "BAD-ROLE",
            event_info="Protein Fusion: in frame; order unavailable",
            connection_type="3to3",
        ),
        _call("NO-BREAKPOINT", breakpoint=None),
    ]

    run = analyze_structural_variant_calls(
        calls,
        "BRAF",
        "edge-study",
        genome_nexus_client=client,
        n_permutations=5,
    )

    assert run.summary["total_fusions"] == 3
    assert run.summary["mapped_fusions"] == 1
    assert run.summary["skipped_fusions"] == 2
    assert [event.Sample_id for event in run.events] == ["VALID"]
    assert any(
        "BAD-ROLE" in warning and "could not determine" in warning
        for warning in run.warnings
    )
    assert any(
        "NO-BREAKPOINT" in warning and "no genomic breakpoint" in warning
        for warning in run.warnings
    )


def test_zero_in_frame_records_emit_outputs_with_unavailable_statistics(
    tmp_path,
    genome_nexus_canonical_transcript_fixture_path,
):
    client = _genome_nexus_client(genome_nexus_canonical_transcript_fixture_path)
    run = analyze_structural_variant_calls(
        [_call("OUT", event_info="Protein Fusion: out of frame  {KIAA1549:BRAF}")],
        "BRAF",
        "edge-study",
        genome_nexus_client=client,
        n_permutations=5,
    )

    assert run.summary["total_fusions"] == 1
    assert run.summary["in_frame_count"] == 0
    assert run.summary["fisher_p_value"] is None
    assert any("no mapped in-frame" in warning for warning in run.warnings)
    paths = write_outputs(run, tmp_path)
    assert "p=unavailable" in paths["markdown"].read_text()
    assert json.loads(paths["json"].read_text())["summary"]["fisher_p_value"] is None


def test_empty_study_result_is_actionable_and_does_not_query_genome_nexus(tmp_path):
    client = MagicMock(spec=GenomeNexusClient)

    run = analyze_structural_variant_calls(
        [],
        "BRAF",
        "typo-study",
        genome_nexus_client=client,
        n_permutations=5,
    )

    assert run.summary["total_fusions"] == 0
    assert run.summary["fisher_p_value"] is None
    assert any("Verify the gene and study ID" in warning for warning in run.warnings)
    client.fetch_canonical_transcript.assert_not_called()
    paths = write_outputs(run, tmp_path)
    assert all(path.exists() for path in paths.values())
    assert "comparison is not possible" in paths["markdown"].read_text()


def test_cbioportal_network_failure_is_wrapped_with_actionable_context(monkeypatch):
    def fail_fetch(*_args, **_kwargs):
        raise requests.Timeout("timed out")

    monkeypatch.setattr(cbioportal_api, "fetch_structural_variants", fail_fetch)

    with pytest.raises(RealBenchmarkNetworkError) as caught:
        run_real_benchmark("BRAF", "some-study")

    message = str(caught.value)
    assert "cBioPortal request failed" in message
    assert "some-study_structural_variants" in message
    assert "Check the study ID" in message


def test_genome_nexus_network_failure_is_wrapped_with_actionable_context():
    client = MagicMock(spec=GenomeNexusClient)
    client.fetch_canonical_transcript.side_effect = requests.ConnectionError("DNS failed")

    with pytest.raises(RealBenchmarkNetworkError) as caught:
        analyze_structural_variant_calls(
            [_call("SAMPLE")],
            "BRAF",
            "some-study",
            genome_nexus_client=client,
            n_permutations=5,
        )

    message = str(caught.value)
    assert "Genome Nexus/domain lookup failed" in message
    assert "Check network access" in message


def test_real_benchmark_click_command_wires_arguments_and_echoes_summary(
    monkeypatch,
    tmp_path,
):
    fake_run = SimpleNamespace(
        gene_symbol="BRAF",
        summary={
            "total_fusions": 7,
            "mapped_fusions": 6,
            "in_frame_count": 5,
            "kinase_retained_count": 4,
            "fisher_p_value": 0.012345,
        },
        warnings=["Skipped one malformed row"],
    )
    run_mock = MagicMock(return_value=fake_run)
    output_paths = {
        "tsv": Path("out.tsv"),
        "json": Path("out.json"),
        "markdown": Path("out.md"),
    }
    write_mock = MagicMock(return_value=output_paths)
    monkeypatch.setattr(cli, "run_real_benchmark", run_mock)
    monkeypatch.setattr(cli, "write_outputs", write_mock)

    result = CliRunner().invoke(
        cli.main,
        [
            "real-benchmark",
            "BRAF",
            "study-id",
            "--output-dir",
            str(tmp_path),
            "--output-stem",
            "custom",
            "--n-permutations",
            "25",
        ],
    )

    assert result.exit_code == 0, result.output
    run_mock.assert_called_once_with("BRAF", "study-id", n_permutations=25)
    write_mock.assert_called_once_with(fake_run, tmp_path, output_stem="custom")
    assert "Analyzed 7 BRAF fusions; mapped=6, in-frame=5" in result.output
    assert "Fisher p=0.012345" in result.output
    assert "Warning: Skipped one malformed row" in result.output
    assert "markdown: out.md" in result.output


def test_real_benchmark_click_command_renders_network_error_without_traceback(monkeypatch):
    monkeypatch.setattr(
        cli,
        "run_real_benchmark",
        MagicMock(side_effect=RealBenchmarkNetworkError("cBioPortal timed out; retry later")),
    )

    result = CliRunner().invoke(cli.main, ["real-benchmark", "BRAF", "study-id"])

    assert result.exit_code == 1
    assert "Error: cBioPortal timed out; retry later" in result.output
    assert "Traceback" not in result.output


def test_real_benchmark_click_command_explains_unknown_gene_without_traceback():
    result = CliRunner().invoke(cli.main, ["real-benchmark", "SOMEGENE", "study-id"])

    assert result.exit_code == 1
    assert "Error: Unknown gene 'SOMEGENE'" in result.output
    assert "cfh list-genes" in result.output
    assert "Traceback" not in result.output
