from cfh.ingestion import clinical_parser


def test_parse_clinical_sample_skips_comment_lines_and_renames_columns(clinical_sample_fixture):
    df = clinical_parser.parse_clinical_sample(clinical_sample_fixture)
    assert list(df.columns) == ["Patient_id", "Sample_id", "Sequencing_panel_id"]
    assert len(df) == 9


def test_parse_clinical_sample_maps_shared_patient_across_samples(clinical_sample_fixture):
    df = clinical_parser.parse_clinical_sample(clinical_sample_fixture)
    shared = df[df["Patient_id"] == "PATIENT-001"]
    assert set(shared["Sample_id"]) == {"SAMPLE-001", "SAMPLE-002"}


def test_parse_clinical_patient_skips_comment_lines(clinical_patient_fixture):
    df = clinical_parser.parse_clinical_patient(clinical_patient_fixture)
    assert "Patient_id" in df.columns
    assert len(df) == 8
