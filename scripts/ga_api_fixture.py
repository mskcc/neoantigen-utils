"""Convert a cBioPortal study fragment into API-shaped JSON for the tab's tests.

Usage: ga_api_fixture.py FRAGMENT_DIR STUDY_ID SAMPLE_ID PATIENT_ID OUT_JSON
"""

import argparse
import csv
import json
from pathlib import Path


def read_meta(path):
    meta = {}
    for line in path.read_text().splitlines():
        key, _, value = line.partition(":")
        meta[key.strip()] = value.strip()
    return meta


def read_table(path):
    with open(path, newline="") as fh:
        return list(csv.reader(fh, delimiter="\t"))


def data_row(profile_id, entity, value, ids):
    return {
        "molecularProfileId": f"{ids['studyId']}_{profile_id}",
        "genericAssayStableId": entity,
        "stableId": entity,
        "sampleId": ids["sampleId"],
        "patientId": ids["patientId"],
        "studyId": ids["studyId"],
        "patientLevel": False,
        "uniqueSampleKey": "",
        "uniquePatientKey": "",
        "value": value,
    }


def profile_rows(meta, frag, ids):
    """Yield (data rows, meta rows) for one GENERIC_ASSAY profile."""
    props = [p for p in meta["generic_entity_meta_properties"].split(",") if p]
    header, *body = read_table(frag / meta["data_filename"])
    col = {name: i for i, name in enumerate(header)}
    for needed in [ids["sampleId"], *props]:
        if needed not in col:
            raise ValueError(f"{meta['data_filename']}: column {needed!r} missing from header")
    sample_col = col[ids["sampleId"]]
    data, metas = [], []
    for row in body:
        entity = row[0]
        metas.append(
            {
                "stableId": entity,
                "entityType": "GENERIC_ASSAY",
                "genericEntityMetaProperties": {p: row[col[p]] for p in props},
            }
        )
        if row[sample_col] != "NA":
            data.append(data_row(meta["stable_id"], entity, row[sample_col], ids))
    return data, metas


def clinical_rows(frag, ids):
    rows = [r for r in read_table(frag / "data_clinical_sample.txt") if not r[0].startswith("#")]
    header, *body = rows
    out, matched = [], False
    for row in body:
        if len(row) != len(header):
            raise ValueError(f"data_clinical_sample.txt: row has {len(row)} cells, header has {len(header)}")
        rec = dict(zip(header, row))
        if rec["SAMPLE_ID"] != ids["sampleId"]:
            continue
        matched = True
        for attr, value in rec.items():
            if attr.startswith("HLA_") and value != "NA":
                out.append(
                    {
                        "clinicalAttributeId": attr,
                        "value": value,
                        "sampleId": ids["sampleId"],
                        "patientId": ids["patientId"],
                        "studyId": ids["studyId"],
                        "patientAttribute": False,
                        "uniqueSampleKey": "",
                        "uniquePatientKey": "",
                    }
                )
    if not matched:
        raise ValueError(f"data_clinical_sample.txt: no row for sample {ids['sampleId']!r}")
    return out


def build_fixture(frag, study_id, sample_id, patient_id):
    frag = Path(frag)
    ids = {"studyId": study_id, "sampleId": sample_id, "patientId": patient_id}
    data, metas, seen = [], [], set()
    for meta_path in sorted(frag.glob("meta_*.txt")):
        meta = read_meta(meta_path)
        if meta.get("genetic_alteration_type") != "GENERIC_ASSAY":
            continue
        d, m = profile_rows(meta, frag, ids)
        data.extend(d)
        for row in m:
            if row["stableId"] not in seen:
                seen.add(row["stableId"])
                metas.append(row)
    if not metas:
        raise ValueError(f"no GENERIC_ASSAY meta files found in {frag}")
    if not data:
        raise ValueError(f"no data rows for sample {sample_id!r} in {frag}")
    return {
        "studyId": study_id,
        "sampleId": sample_id,
        "data": data,
        "meta": metas,
        "clinical": clinical_rows(frag, ids),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("fragment_dir")
    ap.add_argument("study_id")
    ap.add_argument("sample_id")
    ap.add_argument("patient_id")
    ap.add_argument("out_json")
    a = ap.parse_args()
    out = build_fixture(a.fragment_dir, a.study_id, a.sample_id, a.patient_id)
    Path(a.out_json).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out_json).write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
