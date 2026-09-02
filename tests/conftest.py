from pathlib import Path

import pytest

FIXTURES_DIR = Path(__file__).parent / "fixtures"
SV_FIXTURE_DIR = FIXTURES_DIR / "sv"
SV_FIXTURE_FILE = SV_FIXTURE_DIR / "data_sv.txt"
CLINICAL_SAMPLE_FIXTURE = FIXTURES_DIR / "clinical" / "data_clinical_sample.txt"
CLINICAL_PATIENT_FIXTURE = FIXTURES_DIR / "clinical" / "data_clinical_patient.txt"
UNIPROT_FIXTURE = FIXTURES_DIR / "uniprot" / "P15056.json"


@pytest.fixture
def sv_fixture_dir() -> Path:
    return SV_FIXTURE_DIR


@pytest.fixture
def sv_fixture_file() -> Path:
    return SV_FIXTURE_FILE


@pytest.fixture
def clinical_sample_fixture() -> Path:
    return CLINICAL_SAMPLE_FIXTURE


@pytest.fixture
def clinical_patient_fixture() -> Path:
    return CLINICAL_PATIENT_FIXTURE


@pytest.fixture
def uniprot_fixture_path() -> Path:
    return UNIPROT_FIXTURE
