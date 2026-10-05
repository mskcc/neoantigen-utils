import pytest
import importlib.util
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent.parent / "scripts" / "ga_api_fixture.py"


def _load():
    spec = importlib.util.spec_from_file_location("ga_api_fixture", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _write_fragment(d):
    (d / "meta_prof_a.txt").write_text(
        "cancer_study_identifier: s\n"
        "genetic_alteration_type: GENERIC_ASSAY\n"
        "stable_id: prof_a\n"
        "data_filename: data_prof_a.txt\n"
        "generic_entity_meta_properties: NAME,GENE\n"
    )
    (d / "data_prof_a.txt").write_text(
        "ENTITY_STABLE_ID\tNAME\tGENE\tS1\n"
        "ent1\tn1\tTP53\t3\n"
        "ent2\tn2\tKRAS\tNA\n"
    )
    (d / "meta_mut.txt").write_text(
        "genetic_alteration_type: MUTATION_EXTENDED\n"
        "data_filename: data_mut.txt\n"
    )
    (d / "data_mut.txt").write_text("should\tbe\tignored\n")
    (d / "data_clinical_sample.txt").write_text(
        "#a\tb\tc\n#a\tb\tc\n#a\tb\tc\n#a\tb\tc\n"
        "PATIENT_ID\tSAMPLE_ID\tHLA_A_1\tHLA_A_2\tOTHER\n"
        "P1\tS1\tA*30:02\tNA\tx\n"
    )


def test_build_fixture(tmp_path):
    _write_fragment(tmp_path)
    out = _load().build_fixture(tmp_path, "study", "S1", "P1")
    assert out["studyId"] == "study" and out["sampleId"] == "S1"
    assert out["data"] == [
        {
            "molecularProfileId": "study_prof_a",
            "genericAssayStableId": "ent1",
            "stableId": "ent1",
            "sampleId": "S1",
            "patientId": "P1",
            "studyId": "study",
            "patientLevel": False,
            "uniqueSampleKey": "",
            "uniquePatientKey": "",
            "value": "3",
        }
    ]
    assert out["meta"] == [
        {
            "stableId": "ent1",
            "entityType": "GENERIC_ASSAY",
            "genericEntityMetaProperties": {"NAME": "n1", "GENE": "TP53"},
        },
        {
            "stableId": "ent2",
            "entityType": "GENERIC_ASSAY",
            "genericEntityMetaProperties": {"NAME": "n2", "GENE": "KRAS"},
        },
    ]
    assert out["clinical"] == [
        {
            "clinicalAttributeId": "HLA_A_1",
            "value": "A*30:02",
            "sampleId": "S1",
            "patientId": "P1",
            "studyId": "study",
            "patientAttribute": False,
            "uniqueSampleKey": "",
            "uniquePatientKey": "",
        }
    ]


def _build(d, sample="S1"):
    return _load().build_fixture(d, "study", sample, "P1")


def test_zero_generic_assay_metas_raises(tmp_path):
    _write_fragment(tmp_path)
    (tmp_path / "meta_prof_a.txt").unlink()
    with pytest.raises(ValueError, match="no GENERIC_ASSAY"):
        _build(tmp_path)


def test_absent_sample_id_raises(tmp_path):
    _write_fragment(tmp_path)
    with pytest.raises(ValueError, match="'NOPE'"):
        _build(tmp_path, "NOPE")


def test_absent_property_column_raises(tmp_path):
    _write_fragment(tmp_path)
    meta = tmp_path / "meta_prof_a.txt"
    meta.write_text(meta.read_text().replace("NAME,GENE", "NAME,GENE,MISSING"))
    with pytest.raises(ValueError, match="MISSING"):
        _build(tmp_path)


def test_absent_clinical_row_raises(tmp_path):
    _write_fragment(tmp_path)
    clin = tmp_path / "data_clinical_sample.txt"
    clin.write_text(clin.read_text().replace("P1\tS1", "P1\tOTHER"))
    with pytest.raises(ValueError, match="no row for sample"):
        _build(tmp_path)


def test_all_na_sample_raises(tmp_path):
    _write_fragment(tmp_path)
    d = tmp_path / "data_prof_a.txt"
    d.write_text(d.read_text().replace("\t3\n", "\tNA\n"))
    with pytest.raises(ValueError, match="no data rows"):
        _build(tmp_path)
