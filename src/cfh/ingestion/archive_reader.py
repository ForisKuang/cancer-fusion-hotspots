"""Read cBioPortal study data from a ``.tar``/``.tar.gz`` archive or a
plain extracted folder, transparently, so downstream parsing is identical
either way.
"""

from __future__ import annotations

import tarfile
import tempfile
from contextlib import contextmanager
from pathlib import Path
from typing import Iterator

import pandas as pd

from cfh.ingestion.sv_parser import parse_sv_file

SV_FILENAME = "data_sv.txt"


@contextmanager
def open_study_directory(path: str | Path) -> Iterator[Path]:
    """Yield a directory containing the study's files.

    If ``path`` is already a directory it is yielded as-is. If it is a
    ``.tar``/``.tar.gz`` archive it is extracted to a temporary directory
    (cleaned up on exit) and that directory is yielded.
    """
    path = Path(path)
    if path.is_dir():
        yield path
        return

    if tarfile.is_tarfile(path):
        with tempfile.TemporaryDirectory(prefix="cfh_archive_") as tmpdir:
            with tarfile.open(path) as tf:
                tf.extractall(tmpdir, filter="data")
            extracted_root = Path(tmpdir)
            entries = list(extracted_root.iterdir())
            if len(entries) == 1 and entries[0].is_dir():
                yield entries[0]
            else:
                yield extracted_root
        return

    raise ValueError(f"{path} is neither a directory nor a recognized tar archive")


def load_sv_dataframe(path: str | Path) -> pd.DataFrame:
    """Parse ``data_sv.txt`` from a directory or archive at ``path``."""
    with open_study_directory(path) as directory:
        sv_path = directory / SV_FILENAME
        return parse_sv_file(sv_path)
