import copy
import json
import os

import pytest

from neoantigen_utils.cbioportal_study import StudyError, build_study, load_sample

ROOT = {
    "clone_id": 0,
    "X": 1.0,
    "x": 0.0,
    "TMB": 0,
    "neoantigen_load": 0,
    "NA_Mut": 0,
    "F_I": 0.0,
    "F_P": 0,
    "children": [
        {
            "clone_id": 1,
            "X": 1.0,
            "x": 1.0,
            "TMB": 5,
            "neoantigen_load": 1,
            "NA_Mut": 1,
            "F_I": -1.5,
            "F_P": 0,
        }
    ],
}
ANNOTATED = {
    "id": "SAMPLE_1",
    "patient": "PATIENT_1",
    "Effective_N": 168.6,
    "sample_trees": [{"score": -12.5, "topology": ROOT}],
    "mutations": [{"id": "1_100_C_G", "gene": "TP53", "missense": 1}],
    "neoantigens": [
        {
            "id": "n1",
            "mutation_id": "1_100_C_G",
            "HLA_gene_id": "HLA-A02:01",
            "sequence": "ALLAAVLAA",
            "WT_sequence": "ALLAALLAA",
            "mutated_position": 6,
            "Kd": 29.0,
            "KdWT": 26.0,
            "R": 0.0,
            "logC": 1.98,
            "logA": -0.11,
            "quality": 0.5,
        }
    ],
}
TREE = {
    "sample_trees": [
        {
            "score": -12.5,
            "topology": {
                "clone_id": 0,
                "clone_mutations": [],
                "children": [{"clone_id": 1, "clone_mutations": ["1_100_C_G"]}],
            },
        }
    ]
}


def _write_sample(tmp_path):
    annotated = tmp_path / "s1_annotated.json"
    tree = tmp_path / "s1.json"
    annotated.write_text(json.dumps(ANNOTATED))
    tree.write_text(json.dumps(TREE))
    return load_sample("SAMPLE_1", "PATIENT_1", str(annotated), str(tree))


def test_load_sample_returns_all_row_sets(tmp_path):
    sample = _write_sample(tmp_path)
    assert len(sample["nodes"]) == 2
    assert len(sample["neoantigens"]) == 1
    assert sample["mutation_clones"][0]["clone_id"] == 1
    assert sample["summary"]["EFFECTIVE_N"] == 168.6


def test_build_study_writes_clinical_file(tmp_path):
    sample = _write_sample(tmp_path)
    out = tmp_path / "study"
    build_study([sample], "study_1", str(out))

    lines = (out / "data_clinical_sample.txt").read_text().strip().split("\n")
    assert lines[0].startswith("#Patient Identifier")
    assert lines[2].startswith("#STRING\tSTRING")
    assert lines[4].split("\t")[:2] == ["PATIENT_ID", "SAMPLE_ID"]
    assert lines[5].split("\t")[1] == "SAMPLE_1"


def test_build_study_writes_generic_assay_pairs(tmp_path):
    sample = _write_sample(tmp_path)
    out = tmp_path / "study"
    build_study([sample], "study_1", str(out))

    assert os.path.exists(out / "meta_neoantigen_quality.txt")
    data = (out / "data_neoantigen_quality.txt").read_text()
    assert "TP53_1_100_C_G_ALLAAVLAA_A0201" in data
    assert data.strip().split("\n")[0].endswith("SAMPLE_1")


def test_build_study_hides_clone_profiles(tmp_path):
    sample = _write_sample(tmp_path)
    out = tmp_path / "study"
    build_study([sample], "study_1", str(out))

    meta = (out / "meta_clone_parent_t1.txt").read_text()
    assert "show_profile_in_analysis_tab: false" in meta
    data = (out / "data_clone_parent_t1.txt").read_text()
    rows = {line.split("\t")[0]: line.split("\t")[-1] for line in data.strip().split("\n")[1:]}
    assert rows["clone_0"] == "-1"
    assert rows["clone_1"] == "0"


def test_build_study_writes_case_list_and_mutation_columns(tmp_path):
    sample = _write_sample(tmp_path)
    out = tmp_path / "study"
    build_study([sample], "study_1", str(out))

    case_list = (out / "case_lists" / "cases_neoantigen.txt").read_text()
    assert "stable_id: study_1_neoantigen" in case_list
    assert "case_list_ids: SAMPLE_1" in case_list

    columns = (out / "data_neoag_mutation_columns.txt").read_text()
    header = columns.strip().split("\n")[0].split("\t")
    assert header[:3] == ["SAMPLE_ID", "mutation_id", "neoag.clone_id_t1"]
    assert "neoag.ccf" in header
    assert "neoag.clonal" in header


def test_neoantigen_profiles_carry_per_profile_sort_order(tmp_path):
    """Kd is an affinity in nM: lower binds stronger, so it must sort ASC, not DESC."""
    sample = _write_sample(tmp_path)
    out = tmp_path / "study"
    build_study([sample], "study_1", str(out))

    # Anchored to the whole line: a bare "DESC" also matches "DESCRIPTION".
    def sort_order(stable_id):
        lines = (out / "meta_{}.txt".format(stable_id)).read_text().split("\n")
        return next(line for line in lines if line.startswith("value_sort_order:"))

    assert sort_order("neoantigen_kd") == "value_sort_order: ASC"
    assert sort_order("neoantigen_kdwt") == "value_sort_order: ASC"
    assert sort_order("neoantigen_quality") == "value_sort_order: DESC"
    assert sort_order("neoantigen_r") == "value_sort_order: DESC"
    assert sort_order("neoantigen_logc") == "value_sort_order: DESC"
    assert sort_order("neoantigen_loga") == "value_sort_order: DESC"


# A mitochondrial call assigned to no clone that still produces neoantigens. Real
# samples carry these: 3OLTS has one such mutation yielding 2 neoantigens.
UNASSIGNED_MUTATION = "MT_1000_A_T"


def _write_unassigned_sample(tmp_path):
    annotated = copy.deepcopy(ANNOTATED)
    annotated["mutations"].append({"id": UNASSIGNED_MUTATION, "gene": "MT-ND1", "missense": 1})
    for index, quality in enumerate([0.2, 0.9]):
        annotated["neoantigens"].append(
            {
                "id": "mt{}".format(index),
                "mutation_id": UNASSIGNED_MUTATION,
                "HLA_gene_id": "HLA-B07:02",
                "sequence": "MPPLLAAA" + str(index),
                "WT_sequence": "MPPLLAAAA",
                "mutated_position": 3,
                "Kd": 40.0,
                "KdWT": 38.0,
                "R": 0.0,
                "logC": 1.0,
                "logA": -0.2,
                "quality": quality,
            }
        )
    annotated_path = tmp_path / "s2_annotated.json"
    tree_path = tmp_path / "s2.json"
    annotated_path.write_text(json.dumps(annotated))
    # The tree JSON is unchanged: the MT mutation is in no clone_mutations list.
    tree_path.write_text(json.dumps(TREE))
    return load_sample("SAMPLE_1", "PATIENT_1", str(annotated_path), str(tree_path))


def test_mutation_columns_include_unassigned_neoantigenic_mutation(tmp_path):
    """A mutation in no clone still reaches the MAF join, with unknown clonality."""
    sample = _write_unassigned_sample(tmp_path)
    out = tmp_path / "study"
    build_study([sample], "study_1", str(out))

    text = (out / "data_neoag_mutation_columns.txt").read_text().strip().split("\n")
    header = text[0].split("\t")
    rows = {line.split("\t")[1]: dict(zip(header, line.split("\t"))) for line in text[1:]}

    assert set(rows) == {"1_100_C_G", UNASSIGNED_MUTATION}
    unassigned = rows[UNASSIGNED_MUTATION]
    assert unassigned["neoag.clonal"] == "Indeterminate"
    assert unassigned["neoag.clone_id_t1"] == "NA"
    assert unassigned["neoag.ccf"] == "NA"
    assert unassigned["neoag.n_neoantigens"] == "2"
    assert unassigned["neoag.best_neoantigen_quality"] == "0.9"
    # The assigned mutation is unaffected.
    assert rows["1_100_C_G"]["neoag.clonal"] == "Clonal"
    assert rows["1_100_C_G"]["neoag.clone_id_t1"] == "1"


def test_build_study_rejects_clone_id_above_entity_cap(tmp_path):
    """A clone id past the reserved entity rows would vanish from the matrix silently."""
    sample = _write_sample(tmp_path)
    sample["nodes"][1]["clone_id"] = 64

    with pytest.raises(StudyError) as excinfo:
        build_study([sample], "study_1", str(tmp_path / "study"))
    assert "SAMPLE_1" in str(excinfo.value)
    assert "64" in str(excinfo.value)
