# generatemutfasta

Construct mutated peptide FASTAs from a MAF for neoantigen prediction, including
**alternate (non-canonical) transcripts** via Genome Nexus' `Additional_Transcripts`
column.

Extracted from the `mskcc-omics-workflows` neoantigen pipeline so the logic can
be versioned, unit-tested, and released independently of the Nextflow module.

## What it does

- Builds WT/MT peptide windows for each mutation from its HGVSc via Mutalyzer.
- **Multi-transcript mode** (`--multi_transcript`): additionally builds peptides
  for the alternate transcripts Genome Nexus reports in the MAF
  `Additional_Transcripts` column (`-m extended`), emitting `*.altMUT.fa`,
  `*.altWT.fa`, and a `*.transcript_map.tsv` provenance side-car.
- Resilient to real-world MAFs: per-variant errors are logged/counted/skipped
  rather than aborting the sample; full-form HGVSc is re-queried versionlessly
  for offline (no-network) cache resolution; missing protein consequences are
  skipped instead of crashing.

## Install

```bash
pip install .            # runtime: pandas (+ mutalyzer, provided by the container)
pip install .[dev]       # + pytest
```

Installing exposes the `generateMutFasta.py` console command, preserving the exact
CLI the execution container invokes.

## Usage

```bash
generateMutFasta.py --sample_id <id> --output_dir <dir> --maf_file <maf> [--multi_transcript]
```

`mutalyzer` must be importable at runtime (it is in the neoantigen execution
container). It is imported tolerantly, so the package installs and its tests run
without it.

## Development

```bash
pip install -e .[dev]
pytest
```

The `Additional_Transcripts` column format (from Genome Nexus `-m extended`) is:
`Transcript_ID,Hugo_Symbol,HGVSp_Short,HGVSc,Variant_Classification`, semicolon-
separated, canonical transcript excluded (it lives in the main MAF columns).
