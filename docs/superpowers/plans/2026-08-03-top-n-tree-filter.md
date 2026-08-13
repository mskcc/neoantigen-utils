# Top-N Tree Filtering by Log Likelihood Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Only feed the top-N (default 10) best-fitting PhyloWGS trees, ranked by `llh`, into `generate_input.py`'s output instead of every tree in `summ.json`.

**Architecture:** Add one pure helper function, `select_top_trees(trees, n)`, that sorts tree keys by descending `llh` and slices to `n`. Wire it into `main()`'s existing tree loop and add a `--top_n_trees` CLI argument (default `10`). No other function changes.

**Tech Stack:** Python 3.12, pytest.

## Global Constraints

- "Top" trees means highest `llh` value (PhyloWGS convention: `llh` closer to 0 = better fit) — sort descending, not ascending.
- If fewer than N trees exist, return all of them — no error, no padding.
- `--top_n_trees` CLI arg, type `int`, default `10`.
- Follow existing code style in `src/neoantigen_utils/generate_input.py` (plain functions, no classes, google-style-ish docstrings as seen on `compute_purity`).

---

### Task 1: Add and test `select_top_trees`

**Files:**
- Modify: `src/neoantigen_utils/generate_input.py` (add function near `compute_purity`, e.g. after line 50)
- Test: `tests/test_generate_input.py` (add new test class)

**Interfaces:**
- Produces: `select_top_trees(trees: dict, n: int) -> list[str]` — `trees` is a dict keyed by tree-id strings, each value a dict containing `"llh"` (a number). Returns the keys sorted by descending `llh`, truncated to length `min(n, len(trees))`.

- [ ] **Step 1: Write the failing tests**

Add to `tests/test_generate_input.py`:

```python
from neoantigen_utils.generate_input import select_top_trees


class TestSelectTopTrees:
    def test_returns_top_n_sorted_by_descending_llh(self):
        trees = {
            "0": {"llh": -50.0},
            "1": {"llh": -10.0},
            "2": {"llh": -30.0},
            "3": {"llh": -5.0},
        }
        assert select_top_trees(trees, 2) == ["3", "1"]

    def test_returns_all_keys_when_n_exceeds_tree_count(self):
        trees = {"0": {"llh": -2.0}, "1": {"llh": -1.0}}
        assert select_top_trees(trees, 10) == ["1", "0"]

    def test_n_zero_returns_empty_list(self):
        trees = {"0": {"llh": -2.0}, "1": {"llh": -1.0}}
        assert select_top_trees(trees, 0) == []

    def test_empty_trees_returns_empty_list(self):
        assert select_top_trees({}, 10) == []
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_generate_input.py::TestSelectTopTrees -v`
Expected: FAIL with `ImportError: cannot import name 'select_top_trees'`

- [ ] **Step 3: Implement `select_top_trees`**

Add to `src/neoantigen_utils/generate_input.py`, after `compute_purity` (after line 50):

```python
def select_top_trees(trees, n):
    """Return the ``n`` tree keys from ``trees`` with the highest ``llh``.

    PhyloWGS ``llh`` is a log-likelihood; values closer to 0 indicate a
    better-fitting tree, so this sorts descending.

    :param trees: summ.json's ``trees`` dict, keyed by tree id
    :param n: number of trees to keep; if ``n`` exceeds the number of
              available trees, all tree keys are returned
    :return: list of tree keys, sorted by descending llh, length
             ``min(n, len(trees))``
    """
    return sorted(trees, key=lambda t: trees[t]["llh"], reverse=True)[:n]
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `uv run pytest tests/test_generate_input.py::TestSelectTopTrees -v`
Expected: PASS (4 passed)

- [ ] **Step 5: Commit**

```bash
git add src/neoantigen_utils/generate_input.py tests/test_generate_input.py
git commit -m "feat: add select_top_trees helper for llh-based tree filtering"
```

---

### Task 2: Wire filtering into `main()` and add CLI argument

**Files:**
- Modify: `src/neoantigen_utils/generate_input.py:335` (tree loop), `parse_args()` (around line 1317, alongside `--kD_cutoff`)
- Test: `tests/test_generate_input.py` (add integration-style test)

**Interfaces:**
- Consumes: `select_top_trees(trees, n)` from Task 1.
- Produces: `args.top_n_trees` (int, default 10) available on the parsed args used by `main()`.

- [ ] **Step 1: Write the failing test**

This test exercises the loop-selection logic directly (without needing the full `main()` I/O plumbing — MAF/HLA/netMHCpan files aren't relevant to this behavior). Add to `tests/test_generate_input.py`:

```python
class TestTopNTreesArgDefault:
    def test_default_top_n_trees_is_ten(self):
        import sys
        from neoantigen_utils.generate_input import parse_args

        argv = [
            "prog",
            "--maf_file", "x", "--summary_file", "x", "--mutation_file", "x",
            "--tree_directory", "x", "--id", "x", "--patient_id", "x",
            "--cohort", "x", "--HLA_genes", "x",
            "--netMHCpan_MUT_input", "x", "--netMHCpan_WT_input", "x",
        ]
        old_argv = sys.argv
        sys.argv = argv
        try:
            args = parse_args()
        finally:
            sys.argv = old_argv
        assert args.top_n_trees == 10
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest tests/test_generate_input.py::TestTopNTreesArgDefault -v`
Expected: FAIL with `AttributeError: 'Namespace' object has no attribute 'top_n_trees'`

- [ ] **Step 3: Add the CLI argument**

In `src/neoantigen_utils/generate_input.py`, in `parse_args()`, immediately after the existing `--kD_cutoff` argument (around line 1317-1319):

```python
    parser.add_argument(
        "--top_n_trees",
        type=int,
        default=10,
        help="Number of top trees (by llh, highest first) to include, default is 10",
    )
```

- [ ] **Step 4: Run test to verify it passes**

Run: `uv run pytest tests/test_generate_input.py::TestTopNTreesArgDefault -v`
Expected: PASS

- [ ] **Step 5: Wire filtering into the tree loop**

In `src/neoantigen_utils/generate_input.py`, change line 335 from:

```python
    for tree in trees:
```

to:

```python
    for tree in select_top_trees(trees, args.top_n_trees):
```

- [ ] **Step 6: Run the full test file to confirm nothing broke**

Run: `uv run pytest tests/test_generate_input.py -v`
Expected: all tests PASS

- [ ] **Step 7: Commit**

```bash
git add src/neoantigen_utils/generate_input.py tests/test_generate_input.py
git commit -m "feat: filter to top-N trees by llh in generate_input main loop"
```
