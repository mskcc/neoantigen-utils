#!/usr/bin/env python3
"""Export neoantigen pipeline output as a cBioPortal study directory.

Manifest-driven so it runs as a Databricks job over a table of samples, or
locally over a CSV. Portal identifiers come from the manifest, never from
filenames.
"""

import argparse
import csv

from neoantigen_utils.cbioportal_study import build_study, load_sample

VERSION = 1.0


def export_from_manifest(manifest_path, study_id, outdir):
    """Build a study from a CSV of sample_id,patient_id,annotated_json,tree_json."""
    with open(manifest_path, newline="") as handle:
        rows = list(csv.DictReader(handle))

    samples = [load_sample(r["sample_id"], r["patient_id"], r["annotated_json"], r["tree_json"]) for r in rows]
    build_study(samples, study_id, outdir)


def main():
    parser = argparse.ArgumentParser(description="Export neoantigen pipeline output as a cBioPortal study.")
    parser.add_argument("--manifest", required=True, help="CSV: sample_id,patient_id,annotated_json,tree_json")
    parser.add_argument("--study_id", required=True, help="cancer_study_identifier")
    parser.add_argument("--outdir", required=True, help="Directory to write the study into")
    parser.add_argument("-v", "--version", action="version", version="v{}".format(VERSION))
    args = parser.parse_args()

    export_from_manifest(args.manifest, args.study_id, args.outdir)


if __name__ == "__main__":
    main()
