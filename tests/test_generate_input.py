import pytest

from neoantigen_utils.generate_input import compute_purity


def make_tree_data(structure, prevalences):
    """Build a minimal summ.json-shaped tree dict.

    structure:   {node_id_str: [child_id_int, ...]}
    prevalences: {node_id_str: float}  -> cellular_prevalence[0]
    """
    return {
        "structure": structure,
        "populations": {
            nid: {"cellular_prevalence": [p]} for nid, p in prevalences.items()
        },
    }


class TestComputePurity:
    def test_sums_root_children(self):
        tree_data = make_tree_data(
            structure={"0": [1, 2]},
            prevalences={"0": 1.0, "1": 0.5, "2": 0.3},
        )
        assert compute_purity(tree_data) == pytest.approx(0.8)

    def test_single_root_child(self):
        tree_data = make_tree_data(
            structure={"0": [1]},
            prevalences={"0": 1.0, "1": 0.9},
        )
        assert compute_purity(tree_data) == pytest.approx(0.9)

    def test_missing_root_key_raises(self):
        tree_data = make_tree_data(structure={}, prevalences={"0": 1.0})
        with pytest.raises(ValueError):
            compute_purity(tree_data)

    def test_empty_root_children_raises(self):
        tree_data = make_tree_data(structure={"0": []}, prevalences={"0": 1.0})
        with pytest.raises(ValueError):
            compute_purity(tree_data)

    def test_zero_prevalence_children_raises(self):
        tree_data = make_tree_data(
            structure={"0": [1]},
            prevalences={"0": 1.0, "1": 0.0},
        )
        with pytest.raises(ValueError):
            compute_purity(tree_data)
