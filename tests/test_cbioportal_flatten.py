import pytest

from neoantigen_utils.cbioportal_flatten import (
    FlattenError,
    flatten_mutation_clones,
    flatten_neoantigens,
    flatten_tree_nodes,
    flatten_tree_scores,
    validate_tree_nodes,
)


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


def test_flatten_tree_nodes_tolerates_absent_new_x():
    rows = flatten_tree_nodes(_tree(ROOT), "SAMPLE_1")
    assert rows[1]["new_x"] is None


def test_flatten_tree_nodes_points_each_child_at_its_own_parent():
    # Real trees reach four levels and branch; a single-child chain would not
    # distinguish a correct parent pointer from one that reuses the last node.
    grandchild = dict(LEAF, clone_id=3)
    branching = dict(
        ROOT,
        children=[
            dict(LEAF, clone_id=1, children=[grandchild]),
            dict(LEAF, clone_id=2),
        ],
    )
    rows = flatten_tree_nodes(_tree(branching), "SAMPLE_1")
    # The parent map collapses duplicates, so pin the row count too.
    assert len(rows) == 4
    assert {r["clone_id"]: r["parent"] for r in rows} == {0: -1, 1: 0, 3: 1, 2: 0}


def test_flatten_tree_nodes_numbers_trees_from_one():
    data = {
        "sample_trees": [
            {"score": -1.0, "topology": ROOT},
            {"score": -2.0, "topology": ROOT},
        ]
    }
    rows = flatten_tree_nodes(data, "SAMPLE_1")
    assert sorted({r["tree_idx"] for r in rows}) == [1, 2]


ANNOTATED = {
    "sample_trees": [{"score": -12.5, "topology": ROOT}],
    "mutations": [{"id": "1_100_C_G", "gene": "TP53", "missense": 1}],
    "neoantigens": [
        {
            "id": "1_100_C_G_2_9_A0201",
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


def test_flatten_tree_scores_is_one_row_per_tree():
    rows = flatten_tree_scores(ANNOTATED, "SAMPLE_1")
    assert rows == [{"sample_id": "SAMPLE_1", "tree_idx": 1, "loglik": -12.5}]


def test_flatten_neoantigens_joins_gene_from_mutations():
    rows = flatten_neoantigens(ANNOTATED, "SAMPLE_1")
    assert len(rows) == 1
    assert rows[0]["gene"] == "TP53"
    assert rows[0]["quality"] == 0.5
    assert rows[0]["neoantigen_id"] == "1_100_C_G_2_9_A0201"


def test_flatten_neoantigens_tolerates_unknown_mutation():
    data = dict(ANNOTATED, mutations=[])
    assert flatten_neoantigens(data, "SAMPLE_1")[0]["gene"] == ""


def test_flatten_mutation_clones_reads_pre_annotation_tree():
    tree_data = {
        "sample_trees": [
            {
                "score": -1.0,
                "topology": {
                    "clone_id": 0,
                    "clone_mutations": [],
                    "children": [
                        {"clone_id": 1, "clone_mutations": ["1_100_C_G", "2_200_A_T"]},
                    ],
                },
            }
        ]
    }
    rows = flatten_mutation_clones(tree_data, "SAMPLE_1")
    assert sorted(r["mutation_id"] for r in rows) == ["1_100_C_G", "2_200_A_T"]
    assert {r["clone_id"] for r in rows} == {1}
    assert {r["tree_idx"] for r in rows} == {1}


def _deep_topology():
    # Real trees reach depth 6 and carry clone_mutations well below the root.
    return {
        "clone_id": 0,
        "clone_mutations": [],
        "children": [
            {
                "clone_id": 1,
                "clone_mutations": ["1_100_C_G"],
                "children": [
                    {
                        "clone_id": 2,
                        "clone_mutations": ["2_200_A_T"],
                        "children": [
                            {"clone_id": 3, "clone_mutations": ["3_300_G_A"]},
                        ],
                    }
                ],
            }
        ],
    }


def test_flatten_mutation_clones_attributes_deep_nodes():
    rows = flatten_mutation_clones({"sample_trees": [{"score": -1.0, "topology": _deep_topology()}]}, "SAMPLE_1")
    assert {r["mutation_id"]: r["clone_id"] for r in rows} == {
        "1_100_C_G": 1,
        "2_200_A_T": 2,
        "3_300_G_A": 3,
    }
    assert len(rows) == 3


def test_flatten_mutation_clones_numbers_trees_from_one():
    other = {
        "clone_id": 0,
        "clone_mutations": [],
        "children": [{"clone_id": 1, "clone_mutations": ["9_900_T_C"]}],
    }
    tree_data = {
        "sample_trees": [
            {"score": -1.0, "topology": _deep_topology()},
            {"score": -2.0, "topology": other},
        ]
    }
    rows = flatten_mutation_clones(tree_data, "SAMPLE_1")
    assert sorted((r["tree_idx"], r["mutation_id"]) for r in rows) == [
        (1, "1_100_C_G"),
        (1, "2_200_A_T"),
        (1, "3_300_G_A"),
        (2, "9_900_T_C"),
    ]


def _rows(*specs):
    return [{"sample_id": "SAMPLE_1", "tree_idx": 1, "clone_id": c, "parent": p, "x": x} for c, p, x in specs]


def test_validate_accepts_a_well_formed_tree():
    validate_tree_nodes(_rows((0, -1, 0.5), (1, 0, 0.5)))


def test_validate_rejects_two_roots():
    with pytest.raises(FlattenError, match="expected 1 root"):
        validate_tree_nodes(_rows((0, -1, 0.5), (1, -1, 0.5)))


def test_validate_rejects_prevalence_not_summing_to_one():
    with pytest.raises(FlattenError, match="exclusive prevalence"):
        validate_tree_nodes(_rows((0, -1, 0.5), (1, 0, 0.2)))


def test_validate_rejects_duplicate_clone_ids():
    with pytest.raises(FlattenError, match="duplicate clone"):
        validate_tree_nodes(_rows((0, -1, 0.5), (0, 0, 0.5)))


def _tree_rows(tree_idx, *specs):
    return [dict(row, tree_idx=tree_idx) for row in _rows(*specs)]


def test_validate_checks_each_tree_separately():
    # Taken in aggregate these rows sum to 2.0 and carry two roots; the
    # invariants hold per tree, so validation must accept them.
    rows = _tree_rows(1, (0, -1, 0.5), (1, 0, 0.5)) + _tree_rows(2, (0, -1, 0.4), (1, 0, 0.6))
    validate_tree_nodes(rows)


def test_validate_names_the_offending_tree():
    rows = _tree_rows(1, (0, -1, 0.5), (1, 0, 0.5)) + _tree_rows(2, (0, -1, 0.5), (1, -1, 0.5))
    with pytest.raises(FlattenError, match="tree 2"):
        validate_tree_nodes(rows)
