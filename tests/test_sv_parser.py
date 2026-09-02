import numpy as np
import pandas as pd

from cfh.ingestion import sv_parser

EXPECTED_ROW_COUNT = 10


def _input_row_count(path) -> int:
    with open(path) as fh:
        return sum(1 for _ in fh) - 1  # minus header


def test_parses_every_fixture_row_without_raising(sv_fixture_file):
    df = sv_parser.parse_sv_file(sv_fixture_file)
    assert len(df) == _input_row_count(sv_fixture_file) == EXPECTED_ROW_COUNT


def test_source_row_number_is_1_indexed_and_matches_file_order(sv_fixture_file):
    df = sv_parser.parse_sv_file(sv_fixture_file)
    assert list(df["Source_row_number"]) == list(range(1, EXPECTED_ROW_COUNT + 1))
    # Row order is preserved: first row is the explicit-annotation fusion sample.
    assert df.iloc[0]["Sample_Id"] == "SAMPLE-001"


def test_missing_sample_id_row_is_retained_with_warning(sv_fixture_file):
    df = sv_parser.parse_sv_file(sv_fixture_file)
    missing = df[df["Sample_Id"].isna()]
    assert len(missing) == 1
    assert missing.iloc[0]["Parse_warnings"] is not None
    assert "Sample_Id" in missing.iloc[0]["Parse_warnings"]


def test_impprecise_literal_in_numeric_field_does_not_raise(sv_fixture_file):
    df = sv_parser.parse_sv_file(sv_fixture_file)
    row = df[df["Sample_Id"] == "SAMPLE-008"].iloc[0]
    assert row["Tumor_Split_Read_Count"] == "IMPPRECISE"
    assert row["Parse_warnings"] is not None
    assert "Tumor_Split_Read_Count" in row["Parse_warnings"]


def test_missing_annotation_row_has_no_annotation_but_still_parses(sv_fixture_file):
    df = sv_parser.parse_sv_file(sv_fixture_file)
    row = df[df["Sample_Id"] == "SAMPLE-007"].iloc[0]
    assert pd.isna(row["Annotation"])


def test_valid_rows_have_no_warnings(sv_fixture_file):
    df = sv_parser.parse_sv_file(sv_fixture_file)
    row = df[df["Sample_Id"] == "SAMPLE-001"].iloc[0]
    assert pd.isna(row["Parse_warnings"])
    assert isinstance(row["Site1_Position"], (int, np.integer))


def test_output_is_a_dataframe_with_expected_columns(sv_fixture_file):
    df = sv_parser.parse_sv_file(sv_fixture_file)
    assert isinstance(df, pd.DataFrame)
    for column in sv_parser.EXPECTED_COLUMNS:
        assert column in df.columns
