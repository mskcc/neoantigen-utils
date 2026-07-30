# Purity-Normalized Clone Frequencies in `generate_input.py`

**Date:** 2026-07-30
**Status:** Approved, pending implementation plan

## Problem

`generate_input.py` emits a per-clone `X` and `x` for every node of every sample
tree. Today both fields receive the same value — the raw, un-normalized
PhyloWGS cellular prevalence:

```python
newsubtree["X"] = trees[tree]["populations"][str(subTree)]["cellular_prevalence"][0]
newsubtree["x"] = trees[tree]["populations"][str(subTree)]["cellular_prevalence"][0]
newsubtree["new_x"] = 0.0
```

This block appears twice, once in the recursive branch (`generate_input.py:61-67`)
and once in the base case (`generate_input.py:94-100`).

Reference code governing how these values are meant to be derived shows they are
two distinct quantities, neither of which is the raw prevalence.

## Definitions

For a tree with germline root node `0`:

- **Purity** — the sum of the raw cellular prevalences of the root's children:
  `purity = sum(populations[str(r)]["cellular_prevalence"][0] for r in structure["0"])`
- **`X` (inclusive)** — raw prevalence divided by purity. The root clone (`id 0`)
  is pinned to `1.0` rather than computed.
- **`x` (exclusive)** — `X_i - sum(X_c for c in children(i))`. For a leaf, whose
  child set is empty, `x == X`.

`new_x` is a downstream-populated placeholder and remains `0.0`. The reference
code does not touch it.

## Field Mapping

| Our JSON field | Reference quantity |
|---|---|
| `X` | `X` — purity-normalized inclusive prevalence |
| `x` | `Y` — exclusive frequency |
| `new_x` | not set by this logic; stays `0.0` |

## Sample Index

The reference accepts a sample index `ord` for multi-sample reconstructions. This
implementation keeps the existing hardcoded `cellular_prevalence[0]`, matching
current behavior. No new CLI surface is added. Multi-sample support is out of
scope.

## Implementation

### Per-tree purity

Purity is a property of a tree, not a node, so it is computed once per iteration
of the `for tree in trees:` loop (`generate_input.py:292`), before `makeChild` is
called. `makeChild` reads it via closure, as it already does for `trees`, `tree`,
and `treefile`.

```python
populations = trees[tree]["populations"]
roots = trees[tree]["structure"]["0"]
purity = sum(populations[str(r)]["cellular_prevalence"][0] for r in roots)
```

### Single-pass X and x

The reference computes all `X` values, then makes a second pass for `Y`, because
it operates over a flat `nodes` dict. `makeChild` is recursive and populates
children before filling in the parent, so at the point the parent's fields are
assigned, `newsubtree["children"]` already holds fully-populated child dicts.
Both quantities are therefore computed in the existing single traversal; no
second pass is added.

```python
newsubtree["X"] = (
    1.0
    if int(subTree) == 0
    else populations[str(subTree)]["cellular_prevalence"][0] / purity
)
newsubtree["x"] = newsubtree["X"] - sum(c["X"] for c in newsubtree["children"])
newsubtree["new_x"] = 0.0
```

In the base case `children` is empty, so `x == X`, which is what the reference
produces for leaves. One identical snippet replaces both existing blocks rather
than two divergent copies.

## Error Handling

Two guards the reference omits but this code needs:

1. **`purity == 0`** — possible when `structure["0"]` is empty or absent
   (degenerate single-clone tree). Normalizing would raise `ZeroDivisionError`.
   Instead, raise a `ValueError` naming the offending tree.
2. **The existing bare `except`** at `generate_input.py:68-72` and `101-104`
   currently swallows failures and leaves `X`/`x` at `0`. With normalization in
   play, a silent zero is a wrong frequency propagated downstream rather than a
   missing one. The existing `except` structure is left unchanged (surgical
   change); the purity check runs outside `makeChild` as a per-tree
   precondition, so its error is not caught by it.

## Testing

No test currently covers `generate_input.py`. Add `tests/test_generate_input.py`
built on a small hand-authored PhyloWGS-shaped fixture: a three-clone tree with
known prevalences.

Assertions:

- root `X == 1.0`
- each non-root `X` equals `prevalence / purity`
- a parent's `x` equals its `X` minus the sum of its children's `X`
- a leaf's `x` equals its `X`
- `purity == 0` raises `ValueError`

Load-bearing invariant, checked for every node:

```
x_i + sum(X_c for c in children(i)) == X_i
```

## Validation Against Approved Output

Checked against the reference implementation's published data:
`LukszaLab/NeoantigenEditing`, `data/Patient_data/11-LTS/Recurrent/11_LTS_metastasis.json`
— the pre-annotation file, which is the direct analog of this script's output.
(The `_annotated` sibling is downstream and carries additional `TMB`, `F_I`,
`F_P`, `neoantigen_load`, `NA_Mut` fields.)

Across 57 nodes in 5 sample trees:

| Claim | Result |
|---|---|
| Root `X` pinned to `1.0` | Exactly `1.0` in all 5 trees |
| `X` = prevalence / purity | `sum(root children X) == 1.0` exactly |
| `x = X - sum(children X)` | 0 violations across all 57 nodes |
| Leaf `x == X` | Holds |
| Root `x` | `0.0`, which the formula yields as `1.0 - 1.0`; no special case needed |
| `new_x` is a placeholder | `0.0` everywhere, including post-annotation |
| `clone_mutations` id format | `chrom_pos_ref_alt`, matches this script's existing format |

The root's `x` requires no special-casing: the generic formula produces `0.0`
because the root's children's `X` sum to exactly `1.0` after normalization.

## Known Gap: `tilde_x`

The approved input format carries a field this script does not emit:

```
['X', 'children', 'clone_id', 'clone_mutations', 'new_x', 'tilde_x', 'x']
```

`tilde_x` is not a rescale of `x` — both already sum to `1.0` per tree. It is
`x` with an entire subclade zeroed and the remainder renormalized. In tree 0 of
the reference file, clones 7 and 8 have `tilde_x == 0.0` while every surviving
clone shares the ratio `tilde_x / x == 1 / (1 - X_7) == 1.8857`.

Neither the reference snippet governing this work nor any published module in
the `NeoantigenEditing` repository computes `tilde_x`; it is only consumed, in
`predictions_aggregated_loglikelihood_scores.py` (lines 79-105, 178), which
performs the paired primary-versus-recurrent longitudinal analysis. The
single-sample path, `predictions_clones.py:86`, uses plain `x`:

```python
node["predicted_x"] = node["x"] * np.exp(node["fitness"])
```

`tilde_x` therefore appears to originate in an upstream step specific to
paired-sample setups, where the zeroing marks clones absent from the paired
sample. The exclusion rule is unknown and is not inferable from available
sources, so `tilde_x` is deliberately not emitted. This does not block the
single-sample clone-prediction path.

## Out of Scope

- `tilde_x` (see Known Gap above)
- Multi-sample (`ord`) support
- Any change to `new_x` semantics
- Refactoring the existing exception handling in `makeChild`
