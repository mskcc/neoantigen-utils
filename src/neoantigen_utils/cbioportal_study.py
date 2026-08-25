#!/usr/bin/env python3
"""Assemble a loadable cBioPortal study directory from many samples."""

import json
import os

from neoantigen_utils.cbioportal_flatten import (
    flatten_mutation_clones,
    flatten_neoantigens,
    flatten_tree_nodes,
    flatten_tree_scores,
    validate_tree_nodes,
)
from neoantigen_utils.cbioportal_generic_assay import (
    clone_entity_id,
    neoantigen_entity_id,
    render_data,
    render_meta,
)
from neoantigen_utils.cbioportal_study_tables import (
    CLINICAL_ATTRIBUTES,
    CLONE_PROFILES,
    MAX_CLONE_ENTITY,
    NEOANTIGEN_PROFILES,
    N_TREES,
)
from neoantigen_utils.cbioportal_summarize import summarize_sample, truncal_clone


def load_sample(sample_id, patient_id, annotated_path, tree_path):
    """Read one sample's two JSONs and return every row set plus its summary."""
    with open(annotated_path) as handle:
        annotated = json.load(handle)
    with open(tree_path) as handle:
        tree_data = json.load(handle)

    nodes = flatten_tree_nodes(annotated, sample_id)
    validate_tree_nodes(nodes)
    scores = flatten_tree_scores(annotated, sample_id)
    neoantigens = flatten_neoantigens(annotated, sample_id)
    mutation_clones = flatten_mutation_clones(tree_data, sample_id)
    return {
        "sample_id": sample_id,
        "patient_id": patient_id,
        "nodes": nodes,
        "scores": scores,
        "neoantigens": neoantigens,
        "mutation_clones": mutation_clones,
        "summary": summarize_sample(
            nodes,
            scores,
            neoantigens,
            mutation_clones,
            sample_id,
            patient_id,
            annotated.get("Effective_N", "NA"),
        ),
    }


def _write(outdir, filename, text):
    path = os.path.join(outdir, filename)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as handle:
        handle.write(text)


def _render_clinical(samples):
    names = [a[0] for a in CLINICAL_ATTRIBUTES]
    lines = [
        "#" + "\t".join(a[1] for a in CLINICAL_ATTRIBUTES),
        "#" + "\t".join(a[2] for a in CLINICAL_ATTRIBUTES),
        "#" + "\t".join(a[3] for a in CLINICAL_ATTRIBUTES),
        "#" + "\t".join(a[4] for a in CLINICAL_ATTRIBUTES),
        "\t".join(names),
    ]
    for sample in samples:
        lines.append("\t".join(str(sample["summary"].get(n, "NA")) for n in names))
    return "\n".join(lines) + "\n"


def _neoantigen_entities(samples):
    entities = {}
    for sample in samples:
        for row in sample["neoantigens"]:
            entity_id = neoantigen_entity_id(row)
            entities.setdefault(
                entity_id,
                {
                    "ENTITY_STABLE_ID": entity_id,
                    "NAME": row["sequence"],
                    "DESCRIPTION": "{} | {} | {} | WT {} | pos {}".format(
                        row["HLA_gene_id"],
                        row["gene"],
                        row["mutation_id"],
                        row["WT_sequence"],
                        row["mutated_position"],
                    ),
                },
            )
    return [entities[k] for k in sorted(entities)]


def _write_neoantigen_profiles(samples, sample_ids, study_id, outdir):
    entities = _neoantigen_entities(samples)
    for stable_id, field, name, description, sort_order in NEOANTIGEN_PROFILES:
        values = {}
        for sample in samples:
            for row in sample["neoantigens"]:
                values[(neoantigen_entity_id(row), sample["sample_id"])] = row[field]
        data_filename = "data_{}.txt".format(stable_id)
        _write(
            outdir,
            "meta_{}.txt".format(stable_id),
            render_meta(
                study_id,
                stable_id,
                "NEOANTIGEN",
                name,
                description,
                data_filename,
                True,
                value_sort_order=sort_order,
            ),
        )
        _write(outdir, data_filename, render_data(entities, sample_ids, values))


def _write_clone_profiles(samples, sample_ids, study_id, outdir):
    entities = [
        {"ENTITY_STABLE_ID": clone_entity_id(i), "NAME": clone_entity_id(i), "DESCRIPTION": ""}
        for i in range(MAX_CLONE_ENTITY + 1)
    ]
    for tree_idx in range(1, N_TREES + 1):
        for prefix, field, name, description in CLONE_PROFILES:
            stable_id = "{}_t{}".format(prefix, tree_idx)
            values = {}
            for sample in samples:
                for row in sample["nodes"]:
                    if row["tree_idx"] == tree_idx:
                        values[(clone_entity_id(row["clone_id"]), sample["sample_id"])] = row[field]
            if not values:
                continue
            data_filename = "data_{}.txt".format(stable_id)
            _write(
                outdir,
                "meta_{}.txt".format(stable_id),
                render_meta(
                    study_id,
                    stable_id,
                    "CLONE_TREE",
                    "{} (tree {})".format(name, tree_idx),
                    description,
                    data_filename,
                    False,
                ),
            )
            _write(outdir, data_filename, render_data(entities, sample_ids, values))


def _write_case_list(sample_ids, study_id, outdir):
    text = (
        "cancer_study_identifier: {study}\n"
        "stable_id: {study}_neoantigen\n"
        "case_list_name: Samples with neoantigen data\n"
        "case_list_description: Samples processed by the neoantigen pipeline ({n} samples).\n"
        "case_list_ids: {ids}\n"
    ).format(study=study_id, n=len(sample_ids), ids="\t".join(sample_ids))
    _write(outdir, os.path.join("case_lists", "cases_neoantigen.txt"), text)


def _mutation_column_rows(sample):
    """One row per mutation in this sample, over the union of assigned and neoantigenic."""
    top_tree = [r for r in sample["nodes"] if r["tree_idx"] == 1]
    prevalence = {r["clone_id"]: r["X"] for r in top_tree}
    trunk = truncal_clone(top_tree)
    clone_by_tree = {}
    for row in sample["mutation_clones"]:
        clone_by_tree.setdefault(row["mutation_id"], {})[row["tree_idx"]] = row["clone_id"]
    neoantigens_by_mutation = {}
    for row in sample["neoantigens"]:
        neoantigens_by_mutation.setdefault(row["mutation_id"], []).append(row)

    # Union, not just assigned: a mutation in no clone still has neoantigens
    # and still appears in the MAF. Its clonality is Indeterminate.
    rows = []
    for mutation_id in sorted(set(clone_by_tree) | set(neoantigens_by_mutation)):
        per_tree = clone_by_tree.get(mutation_id, {})
        top_clone = per_tree.get(1)
        hits = neoantigens_by_mutation.get(mutation_id, [])
        best = max(hits, key=lambda r: r["quality"]) if hits else None
        row = [sample["sample_id"], mutation_id]
        row += [str(per_tree.get(i, "NA")) for i in range(1, N_TREES + 1)]
        row.append(str(prevalence.get(top_clone, "NA")))
        row.append(_clonality(top_clone, trunk))
        row.append(str(len(hits)))
        row.append(str(best["quality"]) if best else "NA")
        row.append(best["HLA_gene_id"] if best else "NA")
        rows.append(row)
    return rows


def _write_mutation_columns(samples, outdir):
    header = ["SAMPLE_ID", "mutation_id"]
    header += ["neoag.clone_id_t{}".format(i) for i in range(1, N_TREES + 1)]
    header += ["neoag.ccf", "neoag.clonal", "neoag.n_neoantigens", "neoag.best_neoantigen_quality", "neoag.best_hla"]
    lines = ["\t".join(header)]
    for sample in samples:
        lines.extend("\t".join(row) for row in _mutation_column_rows(sample))
    _write(outdir, "data_neoag_mutation_columns.txt", "\n".join(lines) + "\n")


def _clonality(clone_id, trunk):
    """Clonal when the mutation sits on the truncal clone; Indeterminate without a trunk."""
    if trunk is None or clone_id is None:
        return "Indeterminate"
    return "Clonal" if clone_id == trunk["clone_id"] else "Subclonal"


def build_study(samples, study_id, outdir):
    """Write every cBioPortal file for `samples` into `outdir`."""
    os.makedirs(outdir, exist_ok=True)
    sample_ids = [s["sample_id"] for s in samples]

    _write(outdir, "data_clinical_sample.txt", _render_clinical(samples))
    _write(
        outdir,
        "meta_clinical_sample.txt",
        (
            "cancer_study_identifier: {}\n"
            "genetic_alteration_type: CLINICAL\n"
            "datatype: SAMPLE_ATTRIBUTES\n"
            "data_filename: data_clinical_sample.txt\n"
        ).format(study_id),
    )

    _write_neoantigen_profiles(samples, sample_ids, study_id, outdir)
    _write_clone_profiles(samples, sample_ids, study_id, outdir)
    _write_case_list(sample_ids, study_id, outdir)
    _write_mutation_columns(samples, outdir)
