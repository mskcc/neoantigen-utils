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

#### Clone frequency fields

Each node of `sample_trees[].topology` carries two frequencies:

- **`X`** — inclusive clone frequency: the PhyloWGS cellular prevalence divided
  by sample purity. The root clone is pinned to `1.0`, so the root's children
  sum to exactly `1.0`.
- **`x`** — exclusive clone frequency: `X` minus the sum of the children's `X`.
  For a leaf, `x == X`. For the root, `x == 0.0` up to floating-point rounding.

Purity is the sum of the raw cellular prevalences of the root clone's children.
A tree whose purity is `0` is an error, not a tree of zeros.

`new_x` is emitted as `0.0` and populated downstream.

`tilde_x`, present in some NeoantigenEditing inputs, is **not** emitted. It is
only consumed by the paired primary-versus-recurrent longitudinal analysis,
which this package does not target, and the rule determining which clones it
zeroes is not documented in the NeoantigenEditing reference implementation.

The sample index is fixed at `cellular_prevalence[0]`; multi-sample
reconstructions are not supported.

#### Tree selection

`generate_input.py` emits only the best-fitting trees from `summ.json`, ranked
by PhyloWGS `llh` (higher is better) and capped by `--top_n_trees`, which
**defaults to 10**. Prior to VERSION 1.10 every tree in `summ.json` was
emitted; pass a value larger than the tree count to restore that behavior.
Negative values are rejected.

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

`mutalyzer` (HGVS normalization for `generateMutFasta.py`) is an optional,
container-provided dep: it is imported tolerantly / stubbed in tests, so the
package installs and its unit tests run without it.

`pyensembl` (transcript models for `generate_input.py`) is a **required**
dependency. Like `pandas`/`biopython`, it is still imported tolerantly so that
`--help`/`--version` work on a partial install, but it must be present for the
NMD logic to run, and it needs a genome cache at runtime.

## Development

```bash
pip install -e .[dev]
pytest
```

The `Additional_Transcripts` column format (from Genome Nexus `-m extended`) is:
`Transcript_ID,Hugo_Symbol,HGVSp_Short,HGVSc,Variant_Classification`, semicolon-
separated, canonical transcript excluded (it lives in the main MAF columns).
