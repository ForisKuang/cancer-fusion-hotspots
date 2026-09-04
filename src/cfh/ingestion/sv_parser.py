"""Parser for cBioPortal ``data_sv.txt`` structural-variant files.

Never raises on a single malformed row: unparsable values and missing
columns are recorded per-row in ``Parse_warnings`` instead, so the output
row count always equals the input row count.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any

import pandas as pd

EXPECTED_COLUMNS = [
    "Sample_Id",
    "Site1_Hugo_Symbol",
    "Site1_Chromosome",
    "Site1_Position",
    "Site2_Hugo_Symbol",
    "Site2_Chromosome",
    "Site2_Position",
    "Site2_Effect_On_Frame",
    "Tumor_Split_Read_Count",
    "Tumor_Paired_End_Read_Count",
    "SV_Status",
    "NCBI_Build",
    "Connection_Type",
    "Breakpoint_Type",
    "Annotation",
    "Event_Info",
]

_NUMERIC_COLUMNS = {
    "Site1_Position",
    "Site2_Position",
    "Tumor_Split_Read_Count",
    "Tumor_Paired_End_Read_Count",
}

OUTPUT_COLUMNS = EXPECTED_COLUMNS + ["Extra_fields", "Source_row_number", "Parse_warnings"]


def _clean(value: str | None) -> str | None:
    if value is None:
        return None
    value = value.strip()
    return value if value != "" else None


def _coerce_int(value: str) -> int | None:
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _parse_row(raw_row: dict[str, Any], line_number: int) -> dict[str, Any]:
    warnings: list[str] = []
    record: dict[str, Any] = {}

    for column in EXPECTED_COLUMNS:
        value = _clean(raw_row.get(column))
        if value is None:
            # Covers both a header the file never declared and a declared
            # column left blank on this row -- either way the field is
            # missing and that's worth flagging, not just for Sample_Id.
            warnings.append(f"missing {column}")
            record[column] = None
            continue
        if column in _NUMERIC_COLUMNS:
            parsed = _coerce_int(value)
            if parsed is None:
                warnings.append(f"could not parse {column}={value!r} as integer")
                record[column] = value
            else:
                record[column] = parsed
        else:
            record[column] = value

    extra = {
        key: value
        for key, value in raw_row.items()
        if key is not None and key not in EXPECTED_COLUMNS
    }

    # csv.DictReader stashes any tab-delimited values beyond the header's
    # column count under the `None` key as a list, instead of raising or
    # dropping them; keep those too rather than silently discarding them.
    surplus_values = raw_row.get(None)
    if surplus_values:
        extra["_surplus_fields"] = surplus_values
        warnings.append(
            f"row has {len(surplus_values)} surplus field(s) beyond the expected header"
        )

    record["Extra_fields"] = extra or None
    record["Source_row_number"] = line_number
    record["Parse_warnings"] = "; ".join(warnings) if warnings else None
    return record


def parse_sv_file(path: str | Path) -> pd.DataFrame:
    """Parse a tab-delimited cBioPortal ``data_sv.txt`` file.

    ``Source_row_number`` is 1-indexed and matches the order data rows
    appear in the file (the header line is not counted).
    """
    path = Path(path)
    records: list[dict[str, Any]] = []
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for line_number, raw_row in enumerate(reader, start=1):
            try:
                records.append(_parse_row(raw_row, line_number))
            except Exception as exc:  # pragma: no cover - defensive, never drop a row
                records.append(
                    {
                        **{col: None for col in EXPECTED_COLUMNS},
                        "Extra_fields": None,
                        "Source_row_number": line_number,
                        "Parse_warnings": f"failed to parse row: {exc}",
                    }
                )
    return pd.DataFrame.from_records(records, columns=OUTPUT_COLUMNS)
