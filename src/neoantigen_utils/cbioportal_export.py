#!/usr/bin/env python3
"""Export neoantigen pipeline output as a cBioPortal study fragment.

The output is NOT a loadable study directory. It holds the generic-assay
profiles, the clinical sample attributes, a case list and the tidy
intermediates; it does not hold `meta_study.txt` (whose `type_of_cancer` must
be a real OncoTree code chosen by a curator) or `data_clinical_patient.txt`.
Merge it into a study directory that already contains those.

Manifest-driven so it runs as a Databricks job over a table of samples, or
locally over a CSV. Portal identifiers come from the manifest, never from
filenames.
"""

import argparse
import csv

from neoantigen_utils.cbioportal_study import build_study, load_sample

VERSION = 1.0


def export_from_manifest(manifest_path, study_id, outdir):
    """Build a study fragment from a CSV of sample_id,patient_id,annotated_json,tree_json."""
    with open(manifest_path, newline="") as handle:
        rows = list(csv.DictReader(handle))

    samples = [load_sample(r["sample_id"], r["patient_id"], r["annotated_json"], r["tree_json"]) for r in rows]
    build_study(samples, study_id, outdir)


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Export neoantigen pipeline output as a cBioPortal study FRAGMENT: generic-assay "
            "profiles, clinical sample attributes, a case list and tidy intermediates. It is not "
            "loadable on its own -- no meta_study.txt and no data_clinical_patient.txt are written."
        )
    )
    parser.add_argument("--manifest", required=True, help="CSV: sample_id,patient_id,annotated_json,tree_json")
    parser.add_argument("--study_id", required=True, help="cancer_study_identifier")
    parser.add_argument(
        "--outdir",
        required=True,
        help="Directory to write the study fragment into; merge it into a study directory that already has "
        "meta_study.txt and data_clinical_patient.txt",
    )
    parser.add_argument("-v", "--version", action="version", version="v{}".format(VERSION))
    args = parser.parse_args()

    export_from_manifest(args.manifest, args.study_id, args.outdir)


if __name__ == "__main__":
    main()
