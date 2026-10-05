import json

from neoantigen_utils.cbioportal_export import export_from_manifest

# A minimal two-clone sample: root plus one child carrying the only mutation,
# which in turn carries the only neoantigen. Kept local rather than imported
# from test_cbioportal_study so the two modules stay independent.
ANNOTATED = {
    "id": "SAMPLE_1",
    "patient": "PATIENT_1",
    "Effective_N": 168.6,
    "sample_trees": [
        {
            "score": -12.5,
            "topology": {
                "clone_id": 0,
                "X": 1.0,
                "x": 0.0,
                "TMB": 0,
                "neoantigen_load": 0,
                "NA_Mut": 0,
                "F_I": 0.0,
                "F_P": 0,
                "children": [
                    {
                        "clone_id": 1,
                        "X": 1.0,
                        "x": 1.0,
                        "TMB": 5,
                        "neoantigen_load": 1,
                        "NA_Mut": 1,
                        "F_I": -1.5,
                        "F_P": 0,
                    }
                ],
            },
        }
    ],
    "mutations": [{"id": "1_100_C_G", "gene": "TP53", "missense": 1}],
    "HLA_genes": ["A*02:01", "A*03:01"],
    "neoantigens": [
        {
            "id": "n1",
            "mutation_id": "1_100_C_G",
            "HLA_gene_id": "HLA-A02:01",
            "sequence": "ALLAAVLAA",
            "WT_sequence": "ALLAALLAA",
            "mutated_position": 6,
            "Kd": 29.0,
            "KdWT": 26.0,
            "R": 0.0,
            "logC": 1.98,
            "logA": -0.11,
            "quality": 0.5,
        }
    ],
}
TREE = {
    "sample_trees": [
        {
            "score": -12.5,
            "topology": {
                "clone_id": 0,
                "clone_mutations": [],
                "children": [{"clone_id": 1, "clone_mutations": ["1_100_C_G"]}],
            },
        }
    ]
}


def test_export_from_manifest_builds_a_study(tmp_path):
    annotated = tmp_path / "s1_annotated.json"
    tree = tmp_path / "s1.json"
    annotated.write_text(json.dumps(ANNOTATED))
    tree.write_text(json.dumps(TREE))

    manifest = tmp_path / "manifest.csv"
    manifest.write_text(
        "sample_id,patient_id,annotated_json,tree_json\n" "SAMPLE_1,PATIENT_1,{},{}\n".format(annotated, tree)
    )

    out = tmp_path / "study"
    export_from_manifest(str(manifest), "study_1", str(out))

    clinical = (out / "data_clinical_sample.txt").read_text()
    assert "SAMPLE_1" in clinical
    assert (out / "data_neoantigen_quality.txt").exists()
