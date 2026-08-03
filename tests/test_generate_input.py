import pytest

from neoantigen_utils.generate_input import compute_purity, select_top_trees


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


from neoantigen_utils.generate_input import build_topology


# Topology: 0 -> {1 -> {2}, 3}
# Raw prevalences: 1=0.8, 2=0.5, 3=0.2  =>  purity = 0.8 + 0.2 = 1.0
FIXTURE_STRUCTURE = {"0": [1, 3], "1": [2]}
FIXTURE_PREVALENCE = {"0": 1.0, "1": 0.8, "2": 0.5, "3": 0.2}
FIXTURE_MUTATIONS = [
    # (ssm_id, chrom, pos, ref, alt, clone)
    ("s0", "1", 100, "A", "T", "1"),
    ("s1", "2", 200, "C", "G", "2"),
    ("s2", "3", 300, "G", "A", "3"),
]


def make_builder_inputs():
    """Build the five arguments build_topology needs."""
    tree_data = make_tree_data(FIXTURE_STRUCTURE, FIXTURE_PREVALENCE)

    mut_assignments = {nid: {"ssms": []} for nid in ["0", "1", "2", "3"]}
    for ssm, _, _, _, _, clone in FIXTURE_MUTATIONS:
        mut_assignments[clone]["ssms"].append(ssm)
    treefile = {"mut_assignments": mut_assignments}

    mut_data = {
        "ssms": {
            ssm: {"name": f"{c}_{p}_{r}_{a}"}
            for ssm, c, p, r, a, _ in FIXTURE_MUTATIONS
        }
    }
    chrom_pos_dict = {
        f"{c}_{p}_{r}_{a}": {"id": f"{c}_{p}_{r}_{a}"}
        for _, c, p, r, a, _ in FIXTURE_MUTATIONS
    }
    return tree_data, treefile, chrom_pos_dict, mut_data


def walk_nodes(node):
    """Yield every node in a topology dict, parents before children."""
    yield node
    for child in node.get("children", []):
        yield from walk_nodes(child)


class TestBuildTopologyStructure:
    def test_shape_and_mutations_preserved(self):
        tree_data, treefile, chrom_pos_dict, mut_data = make_builder_inputs()

        topology = build_topology(
            "0", True, tree_data, treefile, chrom_pos_dict, mut_data, purity=1.0
        )

        by_id = {n["clone_id"]: n for n in walk_nodes(topology)}
        assert set(by_id) == {0, 1, 2, 3}
        assert [c["clone_id"] for c in by_id[0]["children"]] == [1, 3]
        assert [c["clone_id"] for c in by_id[1]["children"]] == [2]
        assert by_id[2]["children"] == []
        assert by_id[3]["children"] == []

        # root carries no mutations; each other clone carries its own
        assert by_id[0]["clone_mutations"] == []
        assert by_id[1]["clone_mutations"] == ["1_100_A_T"]
        assert by_id[2]["clone_mutations"] == ["2_200_C_G"]
        assert by_id[3]["clone_mutations"] == ["3_300_G_A"]

    def test_new_x_is_zero_everywhere(self):
        tree_data, treefile, chrom_pos_dict, mut_data = make_builder_inputs()

        topology = build_topology(
            "0", True, tree_data, treefile, chrom_pos_dict, mut_data, purity=1.0
        )

        for node in walk_nodes(topology):
            assert node["new_x"] == 0.0

    def test_tilde_x_is_not_emitted(self):
        # Documented known gap: we deliberately do not emit tilde_x.
        tree_data, treefile, chrom_pos_dict, mut_data = make_builder_inputs()

        topology = build_topology(
            "0", True, tree_data, treefile, chrom_pos_dict, mut_data, purity=1.0
        )

        for node in walk_nodes(topology):
            assert "tilde_x" not in node


def assert_frequency_invariants(topology):
    """Check the clone-frequency contract on a built topology.

    - root X is exactly 1.0
    - every node: x == X - sum(children X)
    - leaves: x == X
    """
    assert topology["clone_id"] == 0
    assert topology["X"] == 1.0

    for node in walk_nodes(topology):
        children = node.get("children", [])
        expected_x = node["X"] - sum(c["X"] for c in children)
        assert node["x"] == pytest.approx(expected_x), (
            f"clone {node['clone_id']}: x={node['x']} != X - sum(children X)="
            f"{expected_x}"
        )
        if not children:
            assert node["x"] == pytest.approx(node["X"])


class TestFrequencyNormalization:
    def build(self):
        tree_data, treefile, chrom_pos_dict, mut_data = make_builder_inputs()
        purity = compute_purity(tree_data)
        return build_topology(
            "0", True, tree_data, treefile, chrom_pos_dict, mut_data, purity
        )

    def test_invariants_hold(self):
        assert_frequency_invariants(self.build())

    def test_expected_values(self):
        by_id = {n["clone_id"]: n for n in walk_nodes(self.build())}

        # purity == 0.8 + 0.2 == 1.0, so X equals the raw prevalence here
        assert by_id[0]["X"] == 1.0
        assert by_id[1]["X"] == pytest.approx(0.8)
        assert by_id[2]["X"] == pytest.approx(0.5)
        assert by_id[3]["X"] == pytest.approx(0.2)

        assert by_id[0]["x"] == pytest.approx(0.0)  # 1.0 - (0.8 + 0.2)
        assert by_id[1]["x"] == pytest.approx(0.3)  # 0.8 - 0.5
        assert by_id[2]["x"] == pytest.approx(0.5)  # leaf
        assert by_id[3]["x"] == pytest.approx(0.2)  # leaf

    def test_root_children_X_sum_to_one(self):
        topology = self.build()
        assert sum(c["X"] for c in topology["children"]) == pytest.approx(1.0)

    def test_normalization_divides_by_purity(self):
        # Halve every tumour prevalence: purity halves too, so normalised X is
        # unchanged. This is what distinguishes normalised X from raw prevalence.
        tree_data = make_tree_data(
            FIXTURE_STRUCTURE,
            {"0": 1.0, "1": 0.4, "2": 0.25, "3": 0.1},
        )
        _, treefile, chrom_pos_dict, mut_data = make_builder_inputs()
        purity = compute_purity(tree_data)
        assert purity == pytest.approx(0.5)

        topology = build_topology(
            "0", True, tree_data, treefile, chrom_pos_dict, mut_data, purity
        )
        by_id = {n["clone_id"]: n for n in walk_nodes(topology)}

        assert by_id[1]["X"] == pytest.approx(0.8)
        assert by_id[2]["X"] == pytest.approx(0.5)
        assert by_id[3]["X"] == pytest.approx(0.2)
        assert_frequency_invariants(topology)


class TestSelectTopTrees:
    def test_returns_top_n_sorted_by_descending_llh(self):
        trees = {
            "0": {"llh": -50.0},
            "1": {"llh": -10.0},
            "2": {"llh": -30.0},
            "3": {"llh": -5.0},
        }
        assert select_top_trees(trees, 2) == ["3", "1"]

    def test_returns_all_keys_when_n_exceeds_tree_count(self):
        trees = {"0": {"llh": -2.0}, "1": {"llh": -1.0}}
        assert select_top_trees(trees, 10) == ["1", "0"]

    def test_n_zero_returns_empty_list(self):
        trees = {"0": {"llh": -2.0}, "1": {"llh": -1.0}}
        assert select_top_trees(trees, 0) == []

    def test_empty_trees_returns_empty_list(self):
        assert select_top_trees({}, 10) == []


class TestTopNTreesArgDefault:
    def test_default_top_n_trees_is_ten(self):
        import sys
        from neoantigen_utils.generate_input import parse_args

        argv = [
            "prog",
            "--maf_file", "x", "--summary_file", "x", "--mutation_file", "x",
            "--tree_directory", "x", "--id", "x", "--patient_id", "x",
            "--cohort", "x", "--HLA_genes", "x",
            "--netMHCpan_MUT_input", "x", "--netMHCpan_WT_input", "x",
        ]
        old_argv = sys.argv
        sys.argv = argv
        try:
            args = parse_args()
        finally:
            sys.argv = old_argv
        assert args.top_n_trees == 10
