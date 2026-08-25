#!/usr/bin/env python3
"""Render cBioPortal Generic Assay meta and data files.

Generic Assay is an entity x sample matrix. Clone entities are shared and
bounded; neoantigen entities are keyed on peptide plus HLA so recurrent
neoantigens dedupe across samples.
"""

META_TEMPLATE = """cancer_study_identifier: {study_id}
genetic_alteration_type: GENERIC_ASSAY
generic_assay_type: {assay_type}
datatype: LIMIT-VALUE
stable_id: {stable_id}
profile_name: {profile_name}
profile_description: {profile_description}
data_filename: {data_filename}
show_profile_in_analysis_tab: {show_in_analysis}
value_sort_order: {value_sort_order}
patient_level: {patient_level}
generic_entity_meta_properties: {meta_properties}
"""


def neoantigen_entity_id(row):
    """Stable ID for one neoantigen: gene, mutation, peptide and HLA allele."""
    hla = row["HLA_gene_id"].replace("HLA-", "").replace("*", "").replace(":", "")
    return "{}_{}_{}_{}".format(row["gene"], row["mutation_id"], row["sequence"], hla)


def clone_entity_id(clone_id):
    return "clone_{}".format(clone_id)


def render_meta(
    study_id,
    stable_id,
    assay_type,
    profile_name,
    profile_description,
    data_filename,
    show_in_analysis,
    value_sort_order="DESC",
    patient_level=False,
    meta_properties=("NAME", "DESCRIPTION"),
):
    return META_TEMPLATE.format(
        study_id=study_id,
        assay_type=assay_type,
        stable_id=stable_id,
        profile_name=profile_name,
        profile_description=profile_description,
        data_filename=data_filename,
        show_in_analysis="true" if show_in_analysis else "false",
        value_sort_order=value_sort_order,
        patient_level="true" if patient_level else "false",
        meta_properties=",".join(meta_properties),
    )


def render_data(entities, sample_ids, values, meta_properties=("NAME", "DESCRIPTION")):
    """Pivot `values` into an entity x sample matrix. Missing pairs become NA."""
    lines = ["\t".join(["ENTITY_STABLE_ID"] + list(meta_properties) + list(sample_ids))]
    for entity in entities:
        entity_id = entity["ENTITY_STABLE_ID"]
        row = [entity_id] + [str(entity.get(prop, "")) for prop in meta_properties]
        row += [str(values.get((entity_id, s), "NA")) for s in sample_ids]
        lines.append("\t".join(row))
    return "\n".join(lines) + "\n"
