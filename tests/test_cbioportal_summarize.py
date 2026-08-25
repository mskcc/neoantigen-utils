from neoantigen_utils.cbioportal_summarize import summarize_sample, truncal_clone


def _node(clone_id, parent, x, load, tmb=0, fi=0.0, fp=0, tree_idx=1):
    return {
        "sample_id": "SAMPLE_1",
        "tree_idx": tree_idx,
        "clone_id": clone_id,
        "parent": parent,
        "X": x,
        "x": x,
        "TMB": tmb,
        "neoantigen_load": load,
        "NA_Mut": 0,
        "F_I": fi,
        "F_P": fp,
    }


TRUNK_TREE = [
    _node(0, -1, 0.0, 0),
    _node(1, 0, 0.4, 3, tmb=5, fi=-1.0, fp=1),
    _node(2, 1, 0.6, 7, tmb=9, fi=-2.0, fp=0),
]
BRANCHED_TREE = [
    _node(0, -1, 0.0, 0),
    _node(1, 0, 0.5, 3),
    _node(2, 0, 0.5, 4),
]
SCORES = [{"sample_id": "SAMPLE_1", "tree_idx": 1, "loglik": -2420.51}]
NEOANTIGENS = [{"sample_id": "SAMPLE_1", "mutation_id": "m{}".format(i)} for i in range(10)]
# Every neoantigen's mutation is assigned to a clone.
ASSIGNED = [{"sample_id": "SAMPLE_1", "tree_idx": 1, "mutation_id": "m{}".format(i), "clone_id": 1} for i in range(10)]


def test_truncal_clone_is_the_roots_only_child():
    assert truncal_clone(TRUNK_TREE)["clone_id"] == 1


def test_truncal_clone_is_none_when_root_branches():
    assert truncal_clone(BRANCHED_TREE) is None


def test_summary_weights_by_exclusive_prevalence():
    summary = summarize_sample(TRUNK_TREE, SCORES, NEOANTIGENS, ASSIGNED, "SAMPLE_1", "PATIENT_1", 168.6)
    # 0.4*3 + 0.6*7 = 5.4
    assert summary["CCF_WEIGHTED_NEOANTIGEN_LOAD"] == 5.4
    # 0.4*-1.0 + 0.6*-2.0 = -1.6
    assert summary["CCF_WEIGHTED_FITNESS"] == -1.6


def test_summary_reports_dominant_clone_and_counts():
    summary = summarize_sample(TRUNK_TREE, SCORES, NEOANTIGENS, ASSIGNED, "SAMPLE_1", "PATIENT_1", 168.6)
    assert summary["DOMINANT_CLONE_NEOANTIGEN_LOAD"] == 7
    assert summary["N_CLONES"] == 2
    assert summary["MAX_CLONE_F_P"] == 1
    # Excludes the germline root, whose F_I of 0.0 would otherwise always win.
    assert summary["MAX_CLONE_FITNESS"] == -1.0
    assert summary["TREE_LOGLIK"] == -2420.51
    assert summary["SAMPLE_ID"] == "SAMPLE_1"
    assert summary["PATIENT_ID"] == "PATIENT_1"


def test_summary_splits_truncal_and_subclonal_load():
    summary = summarize_sample(TRUNK_TREE, SCORES, NEOANTIGENS, ASSIGNED, "SAMPLE_1", "PATIENT_1", 168.6)
    assert summary["HAS_TRUNCAL_CLONE"] == "TRUE"
    assert summary["TRUNCAL_NEOANTIGEN_LOAD"] == 3
    assert summary["TOTAL_NEOANTIGEN_LOAD"] == 10
    assert summary["SUBCLONAL_NEOANTIGEN_LOAD"] == 7


def test_summary_counts_unassigned_neoantigens_separately():
    # m9's mutation is in no clone — mitochondrial calls really do this.
    partial = [r for r in ASSIGNED if r["mutation_id"] != "m9"]
    summary = summarize_sample(TRUNK_TREE, SCORES, NEOANTIGENS, partial, "SAMPLE_1", "PATIENT_1", 168.6)
    assert summary["TOTAL_NEOANTIGEN_LOAD"] == 10
    assert summary["UNASSIGNED_NEOANTIGEN_LOAD"] == 1
    assert summary["TRUNCAL_NEOANTIGEN_LOAD"] == 3
    # 10 total - 3 truncal - 1 unknown-clonality, NOT 7.
    assert summary["SUBCLONAL_NEOANTIGEN_LOAD"] == 6


def test_summary_marks_truncal_na_without_a_trunk():
    summary = summarize_sample(BRANCHED_TREE, SCORES, NEOANTIGENS, ASSIGNED, "SAMPLE_1", "PATIENT_1", 1.0)
    assert summary["HAS_TRUNCAL_CLONE"] == "FALSE"
    assert summary["TRUNCAL_NEOANTIGEN_LOAD"] == "NA"
    assert summary["SUBCLONAL_NEOANTIGEN_LOAD"] == "NA"
