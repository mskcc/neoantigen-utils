"""Generic Assay profile writers shared by the cBioPortal study builder."""

import os

from neoantigen_utils.cbioportal_flatten import hla_alleles
from neoantigen_utils.cbioportal_generic_assay import (
    neoantigen_entity_id,
    parse_mutation_id,
    render_data,
    render_meta,
)
from neoantigen_utils.cbioportal_study_tables import (
    MUTATION_CLONE_PROFILE,
    MUTATION_DETECTED_PROFILE,
    MUTATION_META_PROPERTIES,
    NEOANTIGEN_CLONE_PROFILE,
    NEOANTIGEN_META_PROPERTIES,
    N_TREES,
)


class StudyError(ValueError):
    """Raised when a study cannot be emitted without losing data."""


def hla_attributes(alleles):
    """Class I alleles -> HLA_<gene>_<slot> clinical attributes, slots in input order."""
    attrs = {}
    for allele in alleles:
        gene = allele.split("*")[0]
        if gene not in ("A", "B", "C"):
            continue
        if "HLA_{}_2".format(gene) in attrs:
            raise StudyError("more than two HLA-{} alleles".format(gene))
        slot = 2 if "HLA_{}_1".format(gene) in attrs else 1
        attrs["HLA_{}_{}".format(gene, slot)] = allele
    return attrs


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


def _mutation_entities(samples):
    """Union of listed, clone-assigned and neoantigenic mutations; listed gene/missense win."""
    info = {}
    for sample in samples:
        for row in sample["neoantigens"]:
            info.setdefault(row["mutation_id"], {"GENE": row["gene"], "MISSENSE": "NA"})
        for row in sample["mutation_clones"]:
            info.setdefault(row["mutation_id"], {"GENE": "", "MISSENSE": "NA"})
        for row in sample["mutations"]:
            info[row["mutation_id"]] = {"GENE": row["gene"], "MISSENSE": row["missense"]}
    entities = []
    for mutation_id in sorted(info):
        entity = {
            "ENTITY_STABLE_ID": mutation_id,
            "NAME": mutation_id,
            "DESCRIPTION": "{} {}".format(info[mutation_id]["GENE"], mutation_id).strip(),
        }
        entity.update(info[mutation_id])
        entity.update(parse_mutation_id(mutation_id))
        entities.append(entity)
    return entities


def _mutation_values(samples, tree_idx):
    values = {}
    for sample in samples:
        if tree_idx is None:
            values.update({(r["mutation_id"], sample["sample_id"]): 1 for r in sample["mutations"]})
        else:
            values.update(
                {
                    (r["mutation_id"], sample["sample_id"]): r["clone_id"]
                    for r in sample["mutation_clones"]
                    if r["tree_idx"] == tree_idx
                }
            )
    return values


def _write_mutation_profiles(samples, sample_ids, study_id, outdir):
    """Hidden MUTATION profiles: detection, plus one clone assignment per candidate tree."""
    entities = _mutation_entities(samples)
    profiles = [(MUTATION_DETECTED_PROFILE, None)]
    profiles += [(MUTATION_CLONE_PROFILE, t) for t in range(1, N_TREES + 1)]
    for (prefix, name, description), tree_idx in profiles:
        if tree_idx is not None and not any(r["tree_idx"] == tree_idx for s in samples for r in s["nodes"]):
            continue
        stable_id = prefix if tree_idx is None else "{}_t{}".format(prefix, tree_idx)
        label = name if tree_idx is None else "{} (tree {})".format(name, tree_idx)
        data_filename = "data_{}.txt".format(stable_id)
        meta = render_meta(
            study_id,
            stable_id,
            "MUTATION",
            label,
            description,
            data_filename,
            False,
            meta_properties=MUTATION_META_PROPERTIES,
        )
        _write(outdir, "meta_{}.txt".format(stable_id), meta)
        values = _mutation_values(samples, tree_idx)
        _write(outdir, data_filename, render_data(entities, sample_ids, values, MUTATION_META_PROPERTIES))
