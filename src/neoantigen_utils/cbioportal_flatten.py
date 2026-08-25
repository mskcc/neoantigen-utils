#!/usr/bin/env python3
"""Flatten neoantigen pipeline JSON into tidytree-shaped rows.

One row per (sample, tree, clone), carrying a parent pointer rather than nested
children. `parent = -1` marks the root, which is the germline node and holds no
mutations of its own.
"""

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
