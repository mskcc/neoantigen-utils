#!/usr/bin/env python3
"""Flatten neoantigen pipeline JSON into tidytree-shaped rows.

One row per (sample, tree, clone), carrying a parent pointer rather than nested
children. `parent = -1` marks the root, which is the germline node and holds no
mutations of its own.
"""

from collections import defaultdict

# new_x is present on every node; tilde_x only on recurrent samples. Both are
# carried through untouched — the tidy table is the durable artifact and must
# not silently drop upstream fields.
NODE_FIELDS = ("TMB", "neoantigen_load", "NA_Mut", "F_I", "F_P", "new_x", "tilde_x")


def _walk(node, parent, tree_idx, sample_id, rows):
    row = {
        "sample_id": sample_id,
        "tree_idx": tree_idx,
        "clone_id": node["clone_id"],
        "parent": parent,
        "X": node.get("X"),
        "x": node.get("x"),
    }
    for field in NODE_FIELDS:
        row[field] = node.get(field)
    rows.append(row)
    for child in node.get("children", []):
        _walk(child, node["clone_id"], tree_idx, sample_id, rows)


def flatten_tree_nodes(data, sample_id):
    """Return one row per (tree, clone) for every candidate tree in `data`."""
    rows = []
    for tree_idx, tree in enumerate(data["sample_trees"], start=1):
        _walk(tree["topology"], -1, tree_idx, sample_id, rows)
    return rows


NEOANTIGEN_FIELDS = (
    "mutation_id",
    "HLA_gene_id",
    "sequence",
    "WT_sequence",
    "mutated_position",
    "Kd",
    "KdWT",
    "R",
    "logC",
    "logA",
    "quality",
)


def flatten_tree_scores(data, sample_id):
    """Return one row per candidate tree, carrying its log-likelihood."""
    return [
        {"sample_id": sample_id, "tree_idx": tree_idx, "loglik": tree["score"]}
        for tree_idx, tree in enumerate(data["sample_trees"], start=1)
    ]


def flatten_neoantigens(data, sample_id):
    """Return one row per neoantigen, with the gene joined in from `mutations`."""
    genes = {m["id"]: m.get("gene", "") for m in data.get("mutations", [])}
    rows = []
    for neoantigen in data.get("neoantigens", []):
        row = {
            "sample_id": sample_id,
            "neoantigen_id": neoantigen["id"],
            "gene": genes.get(neoantigen.get("mutation_id"), ""),
        }
        for field in NEOANTIGEN_FIELDS:
            row[field] = neoantigen.get(field)
        rows.append(row)
    return rows


def flatten_mutation_clones(tree_data, sample_id):
    """Return one row per (tree, mutation) from the PRE-annotation tree JSON.

    `clone_mutations` exists only in the pre-annotation file; the annotated one
    drops it. It is a duplicate-free subset of the mutation set, not a
    partition — some mutations are assigned to no clone at all.
    """
    rows = []
    for tree_idx, tree in enumerate(tree_data["sample_trees"], start=1):
        stack = [tree["topology"]]
        while stack:
            node = stack.pop()
            for mutation_id in node.get("clone_mutations", []):
                rows.append(
                    {
                        "sample_id": sample_id,
                        "tree_idx": tree_idx,
                        "mutation_id": mutation_id,
                        "clone_id": node["clone_id"],
                    }
                )
            stack.extend(node.get("children", []))
    return rows


class FlattenError(ValueError):
    """Raised when the pipeline output violates an assumed invariant."""


def validate_tree_nodes(rows, tolerance=1e-6):
    """Check the invariants the export relies on. Raises FlattenError."""
    by_tree = defaultdict(list)
    for row in rows:
        by_tree[row["tree_idx"]].append(row)

    for tree_idx, tree_rows in sorted(by_tree.items()):
        roots = [r for r in tree_rows if r["parent"] == -1]
        if len(roots) != 1:
            raise FlattenError("tree {}: expected 1 root, found {}".format(tree_idx, len(roots)))

        clone_ids = [r["clone_id"] for r in tree_rows]
        if len(set(clone_ids)) != len(clone_ids):
            raise FlattenError("tree {}: duplicate clone ids".format(tree_idx))

        total = sum(r["x"] for r in tree_rows)
        if abs(total - 1.0) > tolerance:
            raise FlattenError("tree {}: exclusive prevalence sums to {}, expected 1.0".format(tree_idx, total))
