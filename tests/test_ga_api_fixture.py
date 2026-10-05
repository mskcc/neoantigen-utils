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
