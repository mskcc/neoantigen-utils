# Design: Filter PhyloWGS trees to top N by log-likelihood

## Problem

`generate_input.py::main()` currently iterates over every tree in
`summ_data["trees"]` with no filtering (`generate_input.py:333-353`). PhyloWGS
can emit many candidate trees per sample; only the best-fitting ones (by log
likelihood, `llh`) should feed into the output `sample_trees` list.

## Behavior

- Sort tree keys by `trees[tree]["llh"]` descending (PhyloWGS `llh` is a
  log-likelihood where values closer to 0 indicate a better fit).
- Keep only the top `N` trees (default `N = 10`).
- If fewer than `N` trees exist, keep all of them (no error, no padding).
- Everything downstream (topology building, purity computation, `score`
  field) is unchanged — it simply runs over the filtered set instead of the
  full set.

## CLI

Add `--top_n_trees` (default `10`, type `int`) to `parse_args()`.

## Implementation

Add a standalone function, consistent with existing helpers like
`compute_purity` and `_assign_frequencies`:

```python
def select_top_trees(trees, n):
    """Return the ``n`` tree keys from ``trees`` with the highest ``llh``.

    :param trees: summ.json's ``trees`` dict, keyed by tree id
    :param n: number of trees to keep; if ``n`` exceeds the number of
              available trees, all tree keys are returned
    :return: list of tree keys, sorted by descending llh, length
             ``min(n, len(trees))``
    """
    return sorted(trees, key=lambda t: trees[t]["llh"], reverse=True)[:n]
```

In `main()`, replace:

```python
for tree in trees:
```

with:

```python
for tree in select_top_trees(trees, args.top_n_trees):
```

## Testing

Unit test `select_top_trees` directly in `tests/test_generate_input.py`:

- Given a dict of trees with distinct `llh` values, returns the correct
  top-N keys in descending-llh order.
- Given fewer trees than `n`, returns all keys (still sorted).
- Given `n = 0`, returns an empty list.

No changes needed to existing tests that exercise `main()`/`build_topology`
since they already pass through unfiltered (small) tree sets.
