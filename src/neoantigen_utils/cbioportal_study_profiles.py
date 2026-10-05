"""Generic Assay profile writers shared by the cBioPortal study builder."""

import os

from neoantigen_utils.cbioportal_flatten import hla_alleles
from neoantigen_utils.cbioportal_generic_assay import (
    neoantigen_entity_id,
    render_data,
    render_meta,
)
from neoantigen_utils.cbioportal_study_tables import (
    NEOANTIGEN_CLONE_PROFILE,
    NEOANTIGEN_META_PROPERTIES,
    N_TREES,
)


def _write(outdir, filename, text):
    path = os.path.join(outdir, filename)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as handle:
        handle.write(text)


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
                    "GENE": row["gene"],
                    "MUTATION_ID": row["mutation_id"],
                    "SEQUENCE": row["sequence"],
                    "WT_SEQUENCE": row["WT_sequence"],
                    "HLA": hla_alleles({"HLA_genes": [row["HLA_gene_id"]]})[0],
                    "MUTATED_POSITION": str(row["mutated_position"]),
                },
            )
    return [entities[k] for k in sorted(entities)]


def _write_neoantigen_clone_profiles(samples, sample_ids, study_id, outdir):
    """Hidden per-tree profiles: which clone each neoantigen's mutation sits on."""
    entities = _neoantigen_entities(samples)
    prefix, name, description = NEOANTIGEN_CLONE_PROFILE
    for tree_idx in range(1, N_TREES + 1):
        if not any(row["tree_idx"] == tree_idx for sample in samples for row in sample["nodes"]):
            continue
        values = {}
        for sample in samples:
            clone_of = {
                r["mutation_id"]: r["clone_id"] for r in sample["mutation_clones"] if r["tree_idx"] == tree_idx
            }
            for row in sample["neoantigens"]:
                values[(neoantigen_entity_id(row), sample["sample_id"])] = clone_of.get(row["mutation_id"])
        stable_id = "{}_t{}".format(prefix, tree_idx)
        data_filename = "data_{}.txt".format(stable_id)
        meta = render_meta(
            study_id,
            stable_id,
            "NEOANTIGEN",
            "{} (tree {})".format(name, tree_idx),
            description,
            data_filename,
            False,
            meta_properties=NEOANTIGEN_META_PROPERTIES,
        )
        _write(outdir, "meta_{}.txt".format(stable_id), meta)
        _write(outdir, data_filename, render_data(entities, sample_ids, values, NEOANTIGEN_META_PROPERTIES))
