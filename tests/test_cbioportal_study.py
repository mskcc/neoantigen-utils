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
    "HLA_genes": ["A*02:01", "A*03:01"],
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


def test_load_sample_ranks_trees_by_score_so_tree_1_is_best(tmp_path):
    worse = copy.deepcopy(ANNOTATED["sample_trees"][0])
    worse["score"] = -99.0
    worse["topology"]["children"][0]["TMB"] = 7
    annotated = dict(ANNOTATED, sample_trees=[worse, ANNOTATED["sample_trees"][0]])
    tree = {"sample_trees": [dict(TREE["sample_trees"][0], score=-99.0), TREE["sample_trees"][0]]}
    (tmp_path / "a.json").write_text(json.dumps(annotated))
    (tmp_path / "t.json").write_text(json.dumps(tree))
    sample = load_sample("SAMPLE_1", "PATIENT_1", str(tmp_path / "a.json"), str(tmp_path / "t.json"))
    tree_1 = [r for r in sample["nodes"] if r["tree_idx"] == 1]
    assert [r["TMB"] for r in tree_1] == [0, 5]
    assert sample["scores"][0]["loglik"] == -12.5


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

    # A join input for the MAF rather than a portal file, so it lives under tidy/
    # where a data_*.txt with no meta_*.txt beside it is not a validator finding.
    assert not (out / "data_neoag_mutation_columns.txt").exists()
    columns = (out / "tidy" / "data_neoag_mutation_columns.txt").read_text()
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

    text = (out / "tidy" / "data_neoag_mutation_columns.txt").read_text().strip().split("\n")
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


# SAMPLE_2 shares SAMPLE_1's neoantigen — a recurrent one, which must dedupe into
# a single entity row carrying both samples' values — and adds one of its own,
# which must render NA in SAMPLE_1's column.
SHARED_ENTITY = "TP53_1_100_C_G_ALLAAVLAA_A0201"
UNIQUE_ENTITY = "KRAS_2_200_A_T_KLLPPVVAA_A0201"


def _write_second_sample(tmp_path):
    annotated = copy.deepcopy(ANNOTATED)
    annotated["id"] = "SAMPLE_2"
    annotated["neoantigens"][0]["quality"] = 0.75
    annotated["mutations"].append({"id": "2_200_A_T", "gene": "KRAS", "missense": 1})
    annotated["neoantigens"].append(
        dict(
            annotated["neoantigens"][0],
            id="n2",
            mutation_id="2_200_A_T",
            sequence="KLLPPVVAA",
            quality=0.25,
        )
    )
    annotated_path = tmp_path / "s_two_annotated.json"
    tree_path = tmp_path / "s_two.json"
    annotated_path.write_text(json.dumps(annotated))
    tree_path.write_text(json.dumps(TREE))
    return load_sample("SAMPLE_2", "PATIENT_2", str(annotated_path), str(tree_path))


def test_build_study_over_two_samples(tmp_path):
    """The matrix is a union over samples: shared entities dedupe, absent ones read NA."""
    samples = [_write_sample(tmp_path), _write_second_sample(tmp_path)]
    out = tmp_path / "study"
    build_study(samples, "study_1", str(out))

    data = (out / "data_neoantigen_quality.txt").read_text().rstrip("\n").split("\n")
    header = data[0].split("\t")
    assert header[-2:] == ["SAMPLE_1", "SAMPLE_2"]
    rows = {line.split("\t")[0]: line.split("\t") for line in data[1:]}
    # One row for the recurrent neoantigen, not one per sample.
    assert sorted(rows) == [UNIQUE_ENTITY, SHARED_ENTITY]
    assert rows[SHARED_ENTITY][-2:] == ["0.5", "0.75"]
    # Called in SAMPLE_2 only, so SAMPLE_1's cell is NA rather than absent or empty.
    assert rows[UNIQUE_ENTITY][-2:] == ["NA", "0.25"]

    clinical = (out / "data_clinical_sample.txt").read_text().rstrip("\n").split("\n")
    assert [line.split("\t")[:2] for line in clinical[5:]] == [
        ["PATIENT_1", "SAMPLE_1"],
        ["PATIENT_2", "SAMPLE_2"],
    ]

    # The tidy tables concatenate, keeping every sample's rows.
    nodes = _read_tidy(out, "tree_nodes.tsv")
    assert [line.split("\t")[0] for line in nodes[1:]] == ["SAMPLE_1", "SAMPLE_1", "SAMPLE_2", "SAMPLE_2"]
    assert "SAMPLE_2" in (out / "tidy" / "neoantigens.tsv").read_text()

    assert "case_list_ids: SAMPLE_1\tSAMPLE_2" in (out / "case_lists" / "cases_neoantigen.txt").read_text()


def _read_tidy(out, filename):
    """Split an emitted tidy table without eating trailing empty cells.

    A primary sample's rows end in empty `new_x`/`tilde_x` fields, so `.strip()`
    would take the trailing tabs off the final line along with the newline and
    hide a dropped column.
    """
    return (out / "tidy" / filename).read_text().rstrip("\n").split("\n")


def test_build_study_writes_tidytree_tables(tmp_path):
    sample = _write_sample(tmp_path)
    out = tmp_path / "study"
    build_study([sample], "study_1", str(out))

    nodes = _read_tidy(out, "tree_nodes.tsv")
    assert nodes[0].split("\t") == [
        "sample_id",
        "tree_idx",
        "clone_id",
        "parent",
        "X",
        "x",
        "TMB",
        "neoantigen_load",
        "NA_Mut",
        "F_I",
        "F_P",
        "new_x",
        "tilde_x",
    ]
    assert nodes[1].split("\t")[:4] == ["SAMPLE_1", "1", "0", "-1"]
    # tilde_x is absent on a primary sample: an empty cell, not the string "None".
    assert nodes[1].split("\t")[-2:] == ["", ""]
    # The LAST data row must still carry every column. Trailing empty cells are
    # where a dropped column hides, and the last row is where a parse loses them.
    assert len(nodes[-1].split("\t")) == len(nodes[0].split("\t"))

    scores = _read_tidy(out, "tree_scores.tsv")
    assert scores[0].split("\t") == ["sample_id", "tree_idx", "loglik"]
    assert scores[1].split("\t")[2] == "-12.5"

    neoantigens = (out / "tidy" / "neoantigens.tsv").read_text()
    assert "ALLAAVLAA" in neoantigens

    clones = _read_tidy(out, "mutation_clones.tsv")
    assert clones[0].split("\t") == ["sample_id", "tree_idx", "mutation_id", "clone_id"]


def test_build_study_writes_tree_score_profile(tmp_path):
    sample = _write_sample(tmp_path)
    out = tmp_path / "study"
    build_study([sample], "study_1", str(out))

    meta = (out / "meta_clone_tree_score.txt").read_text()
    assert "value_sort_order: DESC" in meta
    assert "show_profile_in_analysis_tab: false" in meta

    data = (out / "data_clone_tree_score.txt").read_text().strip().split("\n")
    rows = {line.split("\t")[0]: line.split("\t")[-1] for line in data[1:]}
    assert rows["tree_1"] == "-12.5"
    assert rows["tree_2"] == "NA"


def test_build_study_rejects_clone_id_above_entity_cap(tmp_path):
    """A clone id past the reserved entity rows would vanish from the matrix silently."""
    sample = _write_sample(tmp_path)
    sample["nodes"][1]["clone_id"] = 64

    with pytest.raises(StudyError) as excinfo:
        build_study([sample], "study_1", str(tmp_path / "study"))
    assert "SAMPLE_1" in str(excinfo.value)
    assert "64" in str(excinfo.value)


def test_build_study_rejects_negative_clone_id(tmp_path):
    """clone_-1 has no entity row either, so the value vanishes exactly as an over-cap id does."""
    sample = _write_sample(tmp_path)
    sample["nodes"][1]["clone_id"] = -1

    with pytest.raises(StudyError) as excinfo:
        build_study([sample], "study_1", str(tmp_path / "study"))
    assert "SAMPLE_1" in str(excinfo.value)
    assert "-1" in str(excinfo.value)


def test_build_study_rejects_more_candidate_trees_than_the_portal_files_hold(tmp_path):
    """Tree 6 would survive only in the tidy tables, with every portal file silently short."""
    sample = _write_sample(tmp_path)
    for tree_idx in range(2, 7):
        sample["nodes"] += [dict(row, tree_idx=tree_idx) for row in sample["nodes"][:2]]

    with pytest.raises(StudyError) as excinfo:
        build_study([sample], "study_1", str(tmp_path / "study"))
    assert "SAMPLE_1" in str(excinfo.value)
    assert "6 candidate trees" in str(excinfo.value)


def test_build_study_accepts_exactly_n_trees(tmp_path):
    """The guard must reject only past the cap: five trees is the current pipeline's output."""
    sample = _write_sample(tmp_path)
    for tree_idx in range(2, 6):
        sample["nodes"] += [dict(row, tree_idx=tree_idx) for row in sample["nodes"][:2]]

    build_study([sample], "study_1", str(tmp_path / "study"))
    assert (tmp_path / "study" / "data_clone_parent_t5.txt").exists()


def test_build_study_rejects_duplicate_sample_ids(tmp_path):
    """A repeated id emits duplicate matrix columns and clinical rows rather than failing."""
    sample = _write_sample(tmp_path)

    with pytest.raises(StudyError, match="duplicate sample_id SAMPLE_1"):
        build_study([sample, sample], "study_1", str(tmp_path / "study"))


def test_generic_assay_renders_an_absent_upstream_field_as_na(tmp_path):
    """An absent field is None on the row; "None" is neither a number nor a null."""
    annotated = copy.deepcopy(ANNOTATED)
    del annotated["neoantigens"][0]["quality"]
    annotated_path = tmp_path / "s3_annotated.json"
    tree_path = tmp_path / "s3.json"
    annotated_path.write_text(json.dumps(annotated))
    tree_path.write_text(json.dumps(TREE))
    sample = load_sample("SAMPLE_1", "PATIENT_1", str(annotated_path), str(tree_path))

    out = tmp_path / "study"
    build_study([sample], "study_1", str(out))

    data = (out / "data_neoantigen_quality.txt").read_text().strip().split("\n")
    row = next(line for line in data[1:] if line.startswith("TP53_"))
    assert row.split("\t")[-1] == "NA"
    assert "None" not in row


def test_load_sample_carries_mutations_and_hla(tmp_path):
    sample = _write_sample(tmp_path)
    assert [m["mutation_id"] for m in sample["mutations"]] == ["1_100_C_G"]
    assert sample["hla"] == ["A*02:01", "A*03:01"]


def test_build_study_writes_mutations_tidy_table(tmp_path):
    sample = _write_sample(tmp_path)
    build_study([sample], "study_1", str(tmp_path / "out"))
    text = (tmp_path / "out" / "tidy" / "mutations.tsv").read_text()
    assert text.splitlines() == ["sample_id\tmutation_id\tgene\tmissense", "SAMPLE_1\t1_100_C_G\tTP53\t1"]
