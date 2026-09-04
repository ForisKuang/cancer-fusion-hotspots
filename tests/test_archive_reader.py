import tarfile

import pandas.testing as pdt

from cfh.ingestion import archive_reader, sv_parser


def test_tar_gz_archive_matches_extracted_folder(sv_fixture_dir, sv_fixture_file, tmp_path):
    archive_path = tmp_path / "study.tar.gz"
    with tarfile.open(archive_path, "w:gz") as tf:
        for file in sv_fixture_dir.iterdir():
            tf.add(file, arcname=file.name)

    from_archive = archive_reader.load_sv_dataframe(archive_path)
    from_folder = sv_parser.parse_sv_file(sv_fixture_file)

    pdt.assert_frame_equal(from_archive, from_folder)


def test_plain_directory_is_read_directly(sv_fixture_dir, sv_fixture_file):
    from_directory = archive_reader.load_sv_dataframe(sv_fixture_dir)
    from_direct_parse = sv_parser.parse_sv_file(sv_fixture_file)
    pdt.assert_frame_equal(from_directory, from_direct_parse)


def test_nested_top_level_folder_in_archive_is_handled(sv_fixture_dir, sv_fixture_file, tmp_path):
    archive_path = tmp_path / "study_nested.tar.gz"
    with tarfile.open(archive_path, "w:gz") as tf:
        for file in sv_fixture_dir.iterdir():
            tf.add(file, arcname=f"msk_impact_50k_2026/{file.name}")

    from_archive = archive_reader.load_sv_dataframe(archive_path)
    from_direct_parse = sv_parser.parse_sv_file(sv_fixture_file)
    pdt.assert_frame_equal(from_archive, from_direct_parse)
