from neoantigen_utils.cbioportal_generic_assay import (
    clone_entity_id,
    neoantigen_entity_id,
    render_data,
    render_meta,
)


def test_neoantigen_entity_id_sanitizes_the_hla_allele():
    row = {
        "gene": "TP53",
        "mutation_id": "1_100_C_G",
        "sequence": "ALLAAVLAA",
        "HLA_gene_id": "HLA-A02:01",
    }
    assert neoantigen_entity_id(row) == "TP53_1_100_C_G_ALLAAVLAA_A0201"


def test_neoantigen_entity_id_handles_star_notation():
    row = {
        "gene": "TP53",
        "mutation_id": "1_100_C_G",
        "sequence": "ALLAAVLAA",
        "HLA_gene_id": "HLA-A*02:01",
    }
    assert neoantigen_entity_id(row) == "TP53_1_100_C_G_ALLAAVLAA_A0201"


def test_clone_entity_id():
    assert clone_entity_id(3) == "clone_3"


def test_render_meta_emits_required_fields():
    text = render_meta(
        study_id="study_1",
        stable_id="neoantigen_quality",
        assay_type="NEOANTIGEN",
        profile_name="Neoantigen quality",
        profile_description="Quality from NeoantigenEditing.",
        data_filename="data_neoantigen_quality.txt",
        show_in_analysis=True,
    )
    assert "genetic_alteration_type: GENERIC_ASSAY" in text
    assert "generic_assay_type: NEOANTIGEN" in text
    assert "datatype: LIMIT-VALUE" in text
    assert "stable_id: neoantigen_quality" in text
    assert "show_profile_in_analysis_tab: true" in text
    assert "patient_level: false" in text
    assert "generic_entity_meta_properties: NAME,DESCRIPTION" in text


def test_render_meta_can_hide_a_profile():
    text = render_meta(
        study_id="study_1",
        stable_id="clone_parent_t1",
        assay_type="CLONE_TREE",
        profile_name="Clone parent",
        profile_description="Parent vector.",
        data_filename="data_clone_parent_t1.txt",
        show_in_analysis=False,
    )
    assert "show_profile_in_analysis_tab: false" in text


def test_render_data_pivots_entities_by_sample():
    entities = [
        {"ENTITY_STABLE_ID": "E1", "NAME": "peptide1", "DESCRIPTION": "d1"},
        {"ENTITY_STABLE_ID": "E2", "NAME": "peptide2", "DESCRIPTION": "d2"},
    ]
    values = {("E1", "SAMPLE_1"): 0.5, ("E2", "SAMPLE_2"): 1.25}
    text = render_data(entities, ["SAMPLE_1", "SAMPLE_2"], values)
    lines = text.strip().split("\n")
    assert lines[0] == "ENTITY_STABLE_ID\tNAME\tDESCRIPTION\tSAMPLE_1\tSAMPLE_2"
    assert lines[1] == "E1\tpeptide1\td1\t0.5\tNA"
    assert lines[2] == "E2\tpeptide2\td2\tNA\t1.25"
