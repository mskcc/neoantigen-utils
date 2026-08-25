#!/usr/bin/env python3
"""Aggregate tidytree rows into cBioPortal clinical sample attributes.

Everything is computed on the top-scoring tree. Candidate trees for one sample
disagree about clone count and about whether a trunk exists, so attributes
derived from structure are conditioned on tree 1 by construction.
"""

TOP_TREE = 1


def truncal_clone(tree_rows):
    """The MRCA of the tumor clones, or None when the root branches.

    The root is the germline node; a root with two or more children means the
    reconstruction found no shared trunk.
    """
    root = next(r for r in tree_rows if r["parent"] == -1)
    children = [r for r in tree_rows if r["parent"] == root["clone_id"]]
    return children[0] if len(children) == 1 else None


def unassigned_neoantigen_count(neoantigen_rows, mutation_clone_rows):
    """Neoantigens whose mutation belongs to no clone in the top-scoring tree.

    `clone_mutations` is a duplicate-free subset of the mutation set, not a
    partition — mitochondrial calls are routinely unassigned. Such neoantigens
    have unknown clonality: neither truncal nor subclonal.
    """
    assigned = {r["mutation_id"] for r in mutation_clone_rows if r["tree_idx"] == TOP_TREE}
    return sum(1 for n in neoantigen_rows if n["mutation_id"] not in assigned)


def summarize_sample(node_rows, score_rows, neoantigen_rows, mutation_clone_rows, sample_id, patient_id, effective_n):
    """Return one clinical sample record for `sample_id`."""
    tree = [r for r in node_rows if r["tree_idx"] == TOP_TREE]
    tumor = [r for r in tree if r["parent"] != -1]
    trunk = truncal_clone(tree)
    total_load = len(neoantigen_rows)
    unassigned = unassigned_neoantigen_count(neoantigen_rows, mutation_clone_rows)
    # Extremes are over tumor clones only: the root is the germline node and its
    # F_I of 0.0 would otherwise win MAX_CLONE_FITNESS in every sample.
    dominant = max(tumor, key=lambda r: r["x"])
    loglik = next(s["loglik"] for s in score_rows if s["tree_idx"] == TOP_TREE)

    return {
        "PATIENT_ID": patient_id,
        "SAMPLE_ID": sample_id,
        "CCF_WEIGHTED_NEOANTIGEN_LOAD": round(sum(r["x"] * r["neoantigen_load"] for r in tree), 4),
        "CCF_WEIGHTED_TMB": round(sum(r["x"] * r["TMB"] for r in tree), 4),
        "CCF_WEIGHTED_FITNESS": round(sum(r["x"] * r["F_I"] for r in tree), 4),
        "DOMINANT_CLONE_FITNESS": dominant["F_I"],
        "DOMINANT_CLONE_NEOANTIGEN_LOAD": dominant["neoantigen_load"],
        "N_CLONES": len(tumor),
        "HAS_TRUNCAL_CLONE": "TRUE" if trunk else "FALSE",
        "TRUNCAL_NEOANTIGEN_LOAD": trunk["neoantigen_load"] if trunk else "NA",
        "TOTAL_NEOANTIGEN_LOAD": total_load,
        "UNASSIGNED_NEOANTIGEN_LOAD": unassigned,
        "SUBCLONAL_NEOANTIGEN_LOAD": (total_load - trunk["neoantigen_load"] - unassigned) if trunk else "NA",
        "MAX_CLONE_FITNESS": max(r["F_I"] for r in tumor),
        "MAX_CLONE_F_P": max(r["F_P"] for r in tumor),
        "TREE_LOGLIK": loglik,
        "EFFECTIVE_N": effective_n,
    }
