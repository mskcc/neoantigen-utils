#!/usr/bin/env python3
"""Assemble this tool's portion of a cBioPortal study from many samples.

What is emitted is a study fragment, not a loadable study: it has no
`meta_study.txt` and no patient-level clinical file, both of which a curator
supplies. See `build_study`.
"""

import json
import os

from neoantigen_utils.cbioportal_flatten import (
    flatten_mutation_clones,
    flatten_mutations,
    flatten_neoantigens,
    flatten_tree_nodes,
    flatten_tree_scores,
    hla_alleles,
    select_top_trees,
    validate_tree_nodes,
)
from neoantigen_utils.cbioportal_generic_assay import (
    clone_entity_id,
    neoantigen_entity_id,
    render_data,
    render_meta,
)
from neoantigen_utils.cbioportal_study_profiles import (
    _neoantigen_entities,
    _write,
    _write_mutation_profiles,
    _write_neoantigen_clone_profiles,
)
from neoantigen_utils.cbioportal_study_tables import (
    CLINICAL_ATTRIBUTES,
    CLONE_PROFILES,
    MAX_CLONE_ENTITY,
    NEOANTIGEN_META_PROPERTIES,
    NEOANTIGEN_PROFILES,
    N_TREES,
    TIDY_TABLES,
    TREE_SCORE_PROFILE,
)
from neoantigen_utils.cbioportal_summarize import summarize_sample, truncal_clone


def load_sample(sample_id, patient_id, annotated_path, tree_path):
    """Read one sample's two JSONs and return every row set plus its summary."""
    with open(annotated_path) as handle:
        annotated = json.load(handle)
    with open(tree_path) as handle:
        tree_data = json.load(handle)

    annotated, tree_data = select_top_trees(annotated, tree_data, N_TREES)
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
        "mutations": flatten_mutations(annotated, sample_id),
        "hla": hla_alleles(annotated),
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
                meta_properties=NEOANTIGEN_META_PROPERTIES,
            ),
        )
        _write(outdir, data_filename, render_data(entities, sample_ids, values, NEOANTIGEN_META_PROPERTIES))


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


def write_tidy_tables(samples, outdir):
    """Write the tidytree-shaped tables the portal files are a projection of."""
    for filename, key, columns in TIDY_TABLES:
        lines = ["\t".join(columns)]
        for sample in samples:
            for row in sample[key]:
                # An absent field is None on the flattened row; write it as an empty
                # cell so the table reads back as null, not the string "None".
                cells = [row.get(c) for c in columns]
                lines.append("\t".join("" if v is None else str(v) for v in cells))
        _write(outdir, os.path.join("tidy", filename), "\n".join(lines) + "\n")


def _write_tree_score_profile(samples, sample_ids, study_id, outdir):
    """Hidden profile holding each candidate tree's log-likelihood, keyed on tree_idx."""
    entities = [
        {"ENTITY_STABLE_ID": "tree_{}".format(i), "NAME": "tree_{}".format(i), "DESCRIPTION": ""}
        for i in range(1, N_TREES + 1)
    ]
    values = {}
    for sample in samples:
        for row in sample["scores"]:
            values[("tree_{}".format(row["tree_idx"]), sample["sample_id"])] = row["loglik"]

    stable_id, name, description = TREE_SCORE_PROFILE
    data_filename = "data_{}.txt".format(stable_id)
    meta = render_meta(study_id, stable_id, "CLONE_TREE", name, description, data_filename, False, "DESC")
    _write(outdir, "meta_{}.txt".format(stable_id), meta)
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
    # A join input for the MAF, not a portal file: a data_*.txt in the study
    # directory with no meta_*.txt beside it is what the validator flags.
    _write(outdir, os.path.join("tidy", "data_neoag_mutation_columns.txt"), "\n".join(lines) + "\n")


def _clonality(clone_id, trunk):
    """Clonal when the mutation sits on the truncal clone; Indeterminate without a trunk."""
    if trunk is None or clone_id is None:
        return "Indeterminate"
    return "Clonal" if clone_id == trunk["clone_id"] else "Subclonal"


class StudyError(ValueError):
    """Raised when a study cannot be emitted without losing data."""


def _check_clone_ids(samples):
    """Clone entity rows cover 0..MAX_CLONE_ENTITY; an id outside it drops silently.

    `render_data` emits one row per entity, so a value keyed on an entity that has
    no row is written nowhere and reported nowhere. A negative id vanishes the same
    way an over-cap one does. Fail loudly instead.
    """
    for sample in samples:
        for row in sample["nodes"]:
            if not 0 <= row["clone_id"] <= MAX_CLONE_ENTITY:
                raise StudyError(
                    "sample {}: clone_id {} outside the clone entity range 0..{}; raise the cap".format(
                        sample["sample_id"], row["clone_id"], MAX_CLONE_ENTITY
                    )
                )


def _check_tree_counts(samples):
    """Every portal projection stops at N_TREES.

    `load_sample` already keeps only the top N_TREES; a sample built any other way
    with more trees would ship tables whose extra trees no portal file carries.
    """
    for sample in samples:
        n_trees = len({row["tree_idx"] for row in sample["nodes"]})
        if n_trees > N_TREES:
            raise StudyError(
                "sample {}: {} candidate trees exceeds N_TREES {}; raise the cap".format(
                    sample["sample_id"], n_trees, N_TREES
                )
            )


def _check_sample_ids(samples):
    """A repeated sample id emits duplicate matrix columns and clinical rows.

    Portal identifiers come from the manifest rather than the filenames, so a
    duplicate is an operator typo — and it corrupts the load instead of aborting it.
    """
    seen = set()
    for sample in samples:
        sample_id = sample["sample_id"]
        if sample_id in seen:
            raise StudyError("duplicate sample_id {}; each sample must appear once in the manifest".format(sample_id))
        seen.add(sample_id)


def build_study(samples, study_id, outdir):
    """Write this tool's portion of a cBioPortal study for `samples` into `outdir`.

    The output is a study FRAGMENT, not a loadable study: generic-assay profiles,
    clinical sample attributes, a case list and the tidy intermediates. It carries
    no `meta_study.txt` and no `data_clinical_patient.txt`, so it must be merged
    into a study directory that already has them.
    """
    _check_sample_ids(samples)
    _check_clone_ids(samples)
    _check_tree_counts(samples)
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
    _write_neoantigen_clone_profiles(samples, sample_ids, study_id, outdir)
    _write_mutation_profiles(samples, sample_ids, study_id, outdir)
    _write_clone_profiles(samples, sample_ids, study_id, outdir)
    _write_tree_score_profile(samples, sample_ids, study_id, outdir)
    write_tidy_tables(samples, outdir)
    _write_case_list(sample_ids, study_id, outdir)
    _write_mutation_columns(samples, outdir)
