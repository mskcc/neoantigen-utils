from neoantigen_utils.cbioportal_flatten import flatten_tree_nodes


def _tree(topology, score=-1.0):
    return {"sample_trees": [{"score": score, "topology": topology}]}


LEAF = {
    "clone_id": 1,
    "X": 0.5,
    "x": 0.5,
    "TMB": 3,
    "neoantigen_load": 2,
    "NA_Mut": 1,
    "F_I": -1.5,
    "F_P": 0,
}
ROOT = {
    "clone_id": 0,
    "X": 1.0,
    "x": 0.5,
    "TMB": 0,
    "neoantigen_load": 0,
    "NA_Mut": 0,
    "F_I": 0.0,
    "F_P": 0,
    "children": [LEAF],
}


def test_flatten_tree_nodes_records_parent_pointers():
    rows = flatten_tree_nodes(_tree(ROOT), "SAMPLE_1")
    assert [r["clone_id"] for r in rows] == [0, 1]
    assert rows[0]["parent"] == -1
    assert rows[1]["parent"] == 0


def test_flatten_tree_nodes_carries_node_attributes():
    rows = flatten_tree_nodes(_tree(ROOT), "SAMPLE_1")
    leaf = rows[1]
    assert leaf["sample_id"] == "SAMPLE_1"
    assert leaf["tree_idx"] == 1
    assert leaf["F_I"] == -1.5
    assert leaf["neoantigen_load"] == 2


def test_flatten_tree_nodes_carries_optional_upstream_fields():
    node = dict(LEAF, new_x=0.25, tilde_x=0.1)
    rows = flatten_tree_nodes(_tree(dict(ROOT, children=[node])), "SAMPLE_1")
    assert rows[1]["new_x"] == 0.25
    assert rows[1]["tilde_x"] == 0.1


def test_flatten_tree_nodes_tolerates_absent_tilde_x():
    # tilde_x appears only on recurrent samples; primaries must still flatten.
    rows = flatten_tree_nodes(_tree(ROOT), "SAMPLE_1")
    assert rows[1]["tilde_x"] is None


def test_flatten_tree_nodes_numbers_trees_from_one():
    data = {
        "sample_trees": [
            {"score": -1.0, "topology": ROOT},
            {"score": -2.0, "topology": ROOT},
        ]
    }
    rows = flatten_tree_nodes(data, "SAMPLE_1")
    assert sorted({r["tree_idx"] for r in rows}) == [1, 2]
