# neoantigen-utils

Helper scripts for the MSK neoantigen pipeline, extracted from the Nextflow
module so the logic can be versioned, unit-tested, and released independently.

Two tools (each installs as its original console command):

### `generateMutFasta.py`
Construct mutated peptide FASTAs from a MAF, including **alternate (non-canonical)
transcripts** via Genome Nexus' `Additional_Transcripts` column.
- Builds WT/MT peptide windows for each mutation from its HGVSc via Mutalyzer.
- **Multi-transcript mode** (`--multi_transcript`): additionally builds peptides
  for the alternate transcripts Genome Nexus reports in the MAF
  `Additional_Transcripts` column (`-m extended`), emitting `*.altMUT.fa`,
  `*.altWT.fa`, and a `*.transcript_map.tsv` provenance side-car.
- Resilient to real-world MAFs: per-variant errors are logged/counted/skipped
  rather than aborting the sample; full-form HGVSc is re-queried versionlessly
  for offline (no-network) cache resolution; missing protein consequences are
  skipped instead of crashing.

### `generate_input.py`
Builds the neoantigen input JSON (summary/mutation) and NMD annotation from a MAF.
Includes a corrected **NMDetective-B** (Lindeboom et al. 2019) implementation:
spliced/CDS-coordinate, strand-aware PTC localization feeding the four nested
decision-tree rules (last exon / <150 nt start-proximal / >407 nt long exon /
50-nt penultimate rule). Uses pyensembl for transcript models.

## Install

```bash
pip install .            # runtime: pandas (+ mutalyzer, provided by the container)
pip install .[dev]       # + pytest
```

Installing exposes both `generateMutFasta.py` and `generate_input.py` console
commands, preserving the exact CLIs the execution container/modules invoke.

## Usage

```bash
generateMutFasta.py --sample_id <id> --output_dir <dir> --maf_file <maf> [--multi_transcript]
generate_input.py --maf_file <maf> --gtf-file <gtf> --cdna-file <cdna> ... 
```

The heavy runtime deps are container-provided and imported tolerantly / stubbed
in tests, so the package installs and its unit tests run without them:
`mutalyzer` (HGVS normalization for `generateMutFasta.py`) and `pyensembl`
(transcript models for `generate_input.py`; needs a genome cache at runtime).

## Development

```bash
pip install -e .[dev]
pytest
```

The `Additional_Transcripts` column format (from Genome Nexus `-m extended`) is:
`Transcript_ID,Hugo_Symbol,HGVSp_Short,HGVSc,Variant_Classification`, semicolon-
separated, canonical transcript excluded (it lives in the main MAF columns).
