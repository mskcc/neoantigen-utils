import pytest

from neoantigen_utils.generate_input import (
    _assign_frequencies,
    compute_purity,
    convert_polysolver_hla,
    select_top_trees,
)


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

    def test_root_child_missing_from_populations_raises_contextualized_error(self):
        # structure lists clone 1 as a root child, but populations has no
        # entry for it at all -- the tree structure and populations disagree.
        tree_data = make_tree_data(structure={"0": [1]}, prevalences={"0": 1.0})
        with pytest.raises(ValueError, match="missing from populations"):
            compute_purity(tree_data)

    def test_root_child_with_no_cellular_prevalence_entries_raises(self):
        # Clone 1 is present in populations, but its cellular_prevalence list
        # is empty -- the IndexError case, not the KeyError case.
        tree_data = make_tree_data(structure={"0": [1]}, prevalences={"0": 1.0})
        tree_data["populations"]["1"] = {"cellular_prevalence": []}
        with pytest.raises(ValueError, match="missing from populations"):
            compute_purity(tree_data)


class TestAssignFrequenciesMalformedPopulations:
    def test_non_root_clone_missing_from_populations_raises_contextualized_error(self):
        node = {"children": []}
        tree_data = {"populations": {}}
        with pytest.raises(ValueError, match="missing from this tree's populations"):
            _assign_frequencies(node, False, "5", tree_data, purity=1.0)

    def test_non_root_clone_with_no_cellular_prevalence_entries_raises(self):
        node = {"children": []}
        tree_data = {"populations": {"5": {"cellular_prevalence": []}}}
        with pytest.raises(ValueError, match="missing from this tree's populations"):
            _assign_frequencies(node, False, "5", tree_data, purity=1.0)

    def test_root_gets_x_one_from_start_flag_without_consulting_populations(self):
        # Regression for the root-ness-desync fix: root-ness now comes purely
        # from the `start` flag build_topology already threads through, not
        # from re-deriving it via int(sub_tree) == 0. Proven here by giving
        # clone "0" no populations entry at all -- if root-ness were still
        # re-derived from sub_tree, this would be indistinguishable from the
        # missing-clone case above and incorrectly raise.
        node = {"children": []}
        tree_data = {"populations": {}}
        _assign_frequencies(node, True, "0", tree_data, purity=1.0)
        assert node["X"] == 1.0


class TestConvertPolysolverHla:
    """convert_polysolver_hla shares hla_string.parse_polysolver_allele's field
    parsing rather than reimplementing it -- these are a regression guard for
    that sharing, mirroring the equivalent hla_string tests."""

    def test_two_digit_second_field(self):
        assert convert_polysolver_hla("hla_a_02_01_01") == "A*02:01"

    def test_three_digit_second_field_is_not_truncated(self):
        # Not a regression guard against this function's own prior behavior --
        # its previous split('_')[2:4] logic already handled this correctly.
        # It's a guard against the *shared* parser drifting in a way that only
        # this call site's tests would catch.
        assert convert_polysolver_hla("hla_b_18_177") == "B*18:177"

    def test_unparseable_entry_raises_value_error(self):
        # The previous inline logic raised a bare, uncontextualized IndexError
        # here (list index out of range on the split result); the shared
        # parser gives a clear ValueError instead.
        with pytest.raises(ValueError, match="Could not parse"):
            convert_polysolver_hla("garbage")


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

    def test_clone_missing_from_mut_assignments_raises(self):
        """summ.json and the tree JSON disagreeing must not emit an empty clone."""
        tree_data, treefile, chrom_pos_dict, mut_data = make_builder_inputs()
        del treefile["mut_assignments"]["2"]

        with pytest.raises(KeyError, match="missing from the tree file"):
            build_topology(
                "0", True, tree_data, treefile, chrom_pos_dict, mut_data, purity=1.0
            )

    def test_unresolvable_ssm_is_dropped_not_fatal(self):
        """An ssm with no MAF record is skipped; the rest of the clone survives."""
        tree_data, treefile, chrom_pos_dict, mut_data = make_builder_inputs()
        treefile["mut_assignments"]["1"]["ssms"].append("s_unknown")

        topology = build_topology(
            "0", True, tree_data, treefile, chrom_pos_dict, mut_data, purity=1.0
        )

        by_id = {n["clone_id"]: n for n in walk_nodes(topology)}
        assert by_id[1]["clone_mutations"] == ["1_100_A_T"]

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
    """PhyloWGS ``llh`` is a log-likelihood; higher (closer to 0) fits better."""

    def test_returns_top_n_sorted_by_descending_llh(self):
        """The best-fitting ``n`` trees come back, best first."""
        trees = {
            "0": {"llh": -50.0},
            "1": {"llh": -10.0},
            "2": {"llh": -30.0},
            "3": {"llh": -5.0},
        }
        assert select_top_trees(trees, 2) == ["3", "1"]

    def test_returns_all_keys_when_n_exceeds_tree_count(self):
        """Asking for more trees than exist keeps all of them, still sorted."""
        trees = {"0": {"llh": -2.0}, "1": {"llh": -1.0}}
        assert select_top_trees(trees, 10) == ["1", "0"]

    def test_n_zero_returns_empty_list(self):
        """``n = 0`` is a valid request for no trees, not an error."""
        trees = {"0": {"llh": -2.0}, "1": {"llh": -1.0}}
        assert select_top_trees(trees, 0) == []

    def test_empty_trees_returns_empty_list(self):
        """No trees in, no trees out."""
        assert select_top_trees({}, 10) == []

    def test_negative_n_raises(self):
        """A negative slice would return the *worst* trees, so reject it."""
        trees = {"0": {"llh": -2.0}, "1": {"llh": -1.0}}
        with pytest.raises(ValueError, match="must be >= 0"):
            select_top_trees(trees, -1)


class TestTopNTreesArg:
    """``--top_n_trees`` caps the emitted trees and must reject negatives."""

    def test_default_top_n_trees_is_ten(self):
        """Omitting the flag keeps the documented default of 10."""
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

    def test_negative_top_n_trees_is_rejected(self):
        """argparse fails fast rather than letting a negative slice through."""
        import argparse

        from neoantigen_utils.generate_input import _non_negative_int

        with pytest.raises(argparse.ArgumentTypeError):
            _non_negative_int("-1")


class TestMafCommentHandling:
    """TEMPO/Genome Nexus MAFs carry a leading ``#version`` line.

    Skipping only the leading ``#`` lines keeps a literal ``#`` inside a data
    row intact; ``pd.read_csv(comment="#")`` would truncate the row there.
    """

    MAF = (
        "#version 2.4\n"
        "Hugo_Symbol\tChromosome\tStart_Position\tVariant_Classification\n"
        "TP53\t17\t7577121\tMissense_Mutation\n"
        "KRAS\t12\t25398284\tMissense_Mutation\n"
    )

    def _write(self, tmp_path):
        maf = tmp_path / "test.maf"
        maf.write_text(self.MAF)
        return maf

    def test_commented_maf_parses(self, tmp_path):
        """The leading ``#version`` line is consumed, not parsed as a header."""
        from neoantigen_utils.generate_input import read_maf

        df = read_maf(self._write(tmp_path))
        assert list(df.columns) == [
            "Hugo_Symbol",
            "Chromosome",
            "Start_Position",
            "Variant_Classification",
        ]
        assert len(df) == 2
        assert df["Hugo_Symbol"].tolist() == ["TP53", "KRAS"]

    def test_uncommented_maf_still_parses(self, tmp_path):
        """A MAF with no comment line is unaffected."""
        from neoantigen_utils.generate_input import read_maf

        maf = tmp_path / "plain.maf"
        maf.write_text(self.MAF.split("\n", 1)[1])
        df = read_maf(maf)
        assert len(df) == 2
        assert df["Hugo_Symbol"].tolist() == ["TP53", "KRAS"]

    def test_hash_inside_a_data_row_is_preserved(self, tmp_path):
        """A literal ``#`` in an annotation field must not truncate the row."""
        from neoantigen_utils.generate_input import read_maf

        maf = tmp_path / "hash.maf"
        maf.write_text(
            "#version 2.4\n"
            "Hugo_Symbol\tNote\tVariant_Classification\n"
            "TP53\tclone #3\tMissense_Mutation\n"
        )
        df = read_maf(maf)
        assert df["Note"].tolist() == ["clone #3"]
        assert df["Variant_Classification"].tolist() == ["Missense_Mutation"]
