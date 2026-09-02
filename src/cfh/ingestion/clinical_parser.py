"""Parsers for cBioPortal clinical data files.

cBioPortal clinical files are tab-delimited with a handful of leading
``#``-prefixed metadata lines before the real header row, and use
SCREAMING_SNAKE_CASE column names that we normalize to this project's
``Field_name`` convention.
"""

from __future__ import annotations

from io import StringIO
from pathlib import Path

import pandas as pd

_COLUMN_RENAME = {
    "PATIENT_ID": "Patient_id",
    "SAMPLE_ID": "Sample_id",
    "SEQUENCING_PANEL_ID": "Sequencing_panel_id",
}


def _read_clinical_file(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    with path.open() as fh:
        data_lines = [line for line in fh if not line.startswith("#")]
    df = pd.read_csv(
        StringIO("".join(data_lines)),
        sep="\t",
        dtype=str,
        keep_default_na=False,
        na_values=[""],
    )
    return df.rename(columns=_COLUMN_RENAME)


def parse_clinical_sample(path: str | Path) -> pd.DataFrame:
    """Parse a ``data_clinical_sample.txt`` file."""
    return _read_clinical_file(path)


def parse_clinical_patient(path: str | Path) -> pd.DataFrame:
    """Parse a ``data_clinical_patient.txt`` file."""
    return _read_clinical_file(path)
