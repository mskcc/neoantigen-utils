# Purity-Normalized Clone Frequencies Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make `generate_input.py` emit purity-normalized inclusive clone frequencies in `X` and exclusive clone frequencies in `x`, instead of writing the same raw PhyloWGS prevalence to both.

**Architecture:** The recursive tree builder `makeChild` is currently a closure inside `main(args)`, reachable only by running the entire pipeline — which requires pyensembl and netMHCpan inputs. Task 2 lifts it to a module-level `build_topology(...)` with its closed-over state passed explicitly, as a behaviour-preserving refactor. Task 3 then implements the frequency change inside that now-testable function. Purity is computed once per tree and passed in. Because the builder populates a node's children before filling in the node itself, both the inclusive and exclusive values are computed in the existing single traversal — no second pass is added.

**Tech Stack:** Python 3.11+, pytest. No new dependencies.

**Spec:** `docs/superpowers/specs/2026-07-30-clone-frequency-normalization-design.md`

## Global Constraints

- Sample index stays hardcoded at `cellular_prevalence[0]`. No new CLI argument.
- `new_x` remains `0.0`. Its semantics are not changed.
- `tilde_x` is NOT emitted. It is a documented known gap; do not attempt to infer it.
- The root clone (`clone_id == 0`) has `X` pinned to `1.0`, not computed.
- Do not alter the `try`/`except` structure inside the tree builder beyond what Task 2's mechanical extraction requires.
- Existing tests must keep passing: `pytest` from the repo root.

## Deviation From the Spec

The spec called this a surgical change that would not restructure `main`. Task 2
adds an extraction the spec did not anticipate. Reason: `makeChild` closes over
`trees`, `tree`, `treefile`, `chrom_pos_dict`, and `mut_data`, so the only way to
reach it is to run `main` end-to-end — and `main` calls `ensembl_load()` at line
349 and reads two netMHCpan TSVs at lines 350-351 before writing its JSON at line
664. The frequency logic cannot be tested without either a pyensembl install plus
a full neoantigen fixture, or this extraction. The extraction is mechanical and
behaviour-preserving.

**Scope of the extraction:** `build_topology` moves out of `main` but stays in
`src/neoantigen_utils/generate_input.py`. No new module is created and the file
is not split. Testability is the sole motivation; file length is not a
consideration here. The spec's "do not refactor the exception handling"
constraint still holds.

The output path is unchanged: `main` continues to write
`args.patient_id + "_" + args.id + "_input.json"` at line 664.

---

### Task 1: Purity helper with a zero-purity guard

**Files:**
- Modify: `src/neoantigen_utils/generate_input.py` (add above `def main(args):` at line 29)
- Test: `tests/test_generate_input.py` (create)

**Interfaces:**
- Produces: `compute_purity(tree_data: dict) -> float`, where `tree_data` is one
  entry of `summ_data["trees"]` — a dict with `"structure"` and `"populations"`
  keys. Returns the sum of `cellular_prevalence[0]` over the root's children.
  Raises `ValueError` when the result is `0`. Tasks 2 and 3 call this.

- [ ] **Step 1: Write the failing tests**

Create `tests/test_generate_input.py`:

```python
import pytest

from neoantigen_utils.generate_input import compute_purity


def make_tree_data(structure, prevalences):
    """Build a minimal summ.json-shaped tree dict.

    structure:   {node_id_str: [child_id_int, ...]}
    prevalences: {node_id_str: float}  -> cellular_prevalence[0]
    """
    return {
        "structure": structure,
        "populations": {
            nid: {"cellular_prevalence": [p]} for nid, p in prevalences.items()
        },
    }


class TestComputePurity:
    def test_sums_root_children(self):
        tree_data = make_tree_data(
            structure={"0": [1, 2]},
            prevalences={"0": 1.0, "1": 0.5, "2": 0.3},
        )
        assert compute_purity(tree_data) == pytest.approx(0.8)

    def test_single_root_child(self):
        tree_data = make_tree_data(
            structure={"0": [1]},
            prevalences={"0": 1.0, "1": 0.9},
        )
        assert compute_purity(tree_data) == pytest.approx(0.9)

    def test_missing_root_key_raises(self):
        tree_data = make_tree_data(structure={}, prevalences={"0": 1.0})
        with pytest.raises(ValueError):
            compute_purity(tree_data)

    def test_empty_root_children_raises(self):
        tree_data = make_tree_data(structure={"0": []}, prevalences={"0": 1.0})
        with pytest.raises(ValueError):
            compute_purity(tree_data)

    def test_zero_prevalence_children_raises(self):
        tree_data = make_tree_data(
            structure={"0": [1]},
            prevalences={"0": 1.0, "1": 0.0},
        )
        with pytest.raises(ValueError):
            compute_purity(tree_data)
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_generate_input.py -v`
Expected: FAIL — `ImportError: cannot import name 'compute_purity'`

- [ ] **Step 3: Write the implementation**

Add to `src/neoantigen_utils/generate_input.py`, above `def main(args):`:

```python
def compute_purity(tree_data):
    """Sum the raw cellular prevalences of the root clone's children.

    PhyloWGS node 0 is the germline population; the tumour clones hanging off
    it sum to the sample purity. Inclusive clone frequencies are normalised by
    this value so that the root's children sum to exactly 1.0.

    :param tree_data: one entry of summ.json's ``trees`` dict, carrying
                      ``structure`` and ``populations``
    :return: purity, strictly greater than 0
    :raises ValueError: if the tree has no root children or they sum to 0
    """
    roots = tree_data["structure"].get("0", [])
    purity = sum(
        tree_data["populations"][str(r)]["cellular_prevalence"][0] for r in roots
    )
    if purity == 0:
        raise ValueError(
            "Cannot normalise clone frequencies: purity is 0 "
            "(no root children, or their cellular prevalences sum to 0)."
        )
    return purity
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_generate_input.py -v`
Expected: PASS, 5 tests

- [ ] **Step 5: Commit**

```bash
git add tests/test_generate_input.py src/neoantigen_utils/generate_input.py
git commit -m "feat: add compute_purity helper for clone frequency normalization"
```

---

### Task 2: Extract `makeChild` to a module-level `build_topology`

Behaviour-preserving refactor. **No frequency logic changes in this task** — `X`
and `x` still both receive the raw prevalence when this task ends. That happens
in Task 3.

**Files:**
- Modify: `src/neoantigen_utils/generate_input.py:31-105` (remove the nested `makeChild`), `:292-303` (update the call site), and add `build_topology` at module level above `main`
- Test: `tests/test_generate_input.py` (extend)

**Interfaces:**
- Produces: `build_topology(sub_tree, start, tree_data, treefile, chrom_pos_dict, mut_data, purity=None) -> dict`
  - `sub_tree`: clone id (int or str) for this node
  - `start`: `True` only for the outermost call, which forces the node id to `0`
  - `tree_data`: one entry of `summ_data["trees"]` (`structure`, `populations`)
  - `treefile`: the per-tree JSON with `mut_assignments`
  - `chrom_pos_dict`: mutation-name to record mapping built earlier in `main`
  - `mut_data`: the parsed mutation JSON with `ssms`
  - `purity`: accepted now, unused until Task 3
  - Returns the topology dict: `clone_id`, `clone_mutations`, `children`, `X`, `x`, `new_x`
  - Task 3 changes the `X`/`x` computation inside this function.

- [ ] **Step 1: Write the characterization test**

Append to `tests/test_generate_input.py`. This pins current behaviour so the
refactor is provably behaviour-preserving:

```python
from neoantigen_utils.generate_input import build_topology


# Topology: 0 -> {1 -> {2}, 3}
# Raw prevalences: 1=0.8, 2=0.5, 3=0.2  =>  purity = 0.8 + 0.2 = 1.0
FIXTURE_STRUCTURE = {"0": [1, 3], "1": [2]}
FIXTURE_PREVALENCE = {"0": 1.0, "1": 0.8, "2": 0.5, "3": 0.2}
FIXTURE_MUTATIONS = [
    # (ssm_id, chrom, pos, ref, alt, clone)
    ("s0", "1", 100, "A", "T", "1"),
    ("s1", "2", 200, "C", "G", "2"),
    ("s2", "3", 300, "G", "A", "3"),
]


def make_builder_inputs():
    """Build the five arguments build_topology needs."""
    tree_data = make_tree_data(FIXTURE_STRUCTURE, FIXTURE_PREVALENCE)

    mut_assignments = {nid: {"ssms": []} for nid in ["0", "1", "2", "3"]}
    for ssm, _, _, _, _, clone in FIXTURE_MUTATIONS:
        mut_assignments[clone]["ssms"].append(ssm)
    treefile = {"mut_assignments": mut_assignments}

    mut_data = {
        "ssms": {
            ssm: {"name": f"{c}_{p}_{r}_{a}"}
            for ssm, c, p, r, a, _ in FIXTURE_MUTATIONS
        }
    }
    chrom_pos_dict = {
        f"{c}_{p}_{r}_{a}": {"id": f"{c}_{p}_{r}_{a}"}
        for _, c, p, r, a, _ in FIXTURE_MUTATIONS
    }
    return tree_data, treefile, chrom_pos_dict, mut_data


def walk_nodes(node):
    """Yield every node in a topology dict, parents before children."""
    yield node
    for child in node.get("children", []):
        yield from walk_nodes(child)


class TestBuildTopologyStructure:
    def test_shape_and_mutations_preserved(self):
        tree_data, treefile, chrom_pos_dict, mut_data = make_builder_inputs()

        topology = build_topology(
            "0", True, tree_data, treefile, chrom_pos_dict, mut_data, purity=1.0
        )

        by_id = {n["clone_id"]: n for n in walk_nodes(topology)}
        assert set(by_id) == {0, 1, 2, 3}
        assert [c["clone_id"] for c in by_id[0]["children"]] == [1, 3]
        assert [c["clone_id"] for c in by_id[1]["children"]] == [2]
        assert by_id[2]["children"] == []
        assert by_id[3]["children"] == []

        # root carries no mutations; each other clone carries its own
        assert by_id[0]["clone_mutations"] == []
        assert by_id[1]["clone_mutations"] == ["1_100_A_T"]
        assert by_id[2]["clone_mutations"] == ["2_200_C_G"]
        assert by_id[3]["clone_mutations"] == ["3_300_G_A"]

    def test_new_x_is_zero_everywhere(self):
        tree_data, treefile, chrom_pos_dict, mut_data = make_builder_inputs()

        topology = build_topology(
            "0", True, tree_data, treefile, chrom_pos_dict, mut_data, purity=1.0
        )

        for node in walk_nodes(topology):
            assert node["new_x"] == 0.0

    def test_tilde_x_is_not_emitted(self):
        # Documented known gap: we deliberately do not emit tilde_x.
        tree_data, treefile, chrom_pos_dict, mut_data = make_builder_inputs()

        topology = build_topology(
            "0", True, tree_data, treefile, chrom_pos_dict, mut_data, purity=1.0
        )

        for node in walk_nodes(topology):
            assert "tilde_x" not in node
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_generate_input.py -v`
Expected: FAIL — `ImportError: cannot import name 'build_topology'`

- [ ] **Step 3: Move the function to module level**

Cut the entire nested `def makeChild(subTree, start):` block from inside `main`
(`generate_input.py:31-105`) and paste it at module level, below
`compute_purity`. Rename to `build_topology`, add the parameters, and change the
closed-over names to the new parameters. The only substantive edits are the
signature, the two recursive calls, and swapping `trees[tree]` for `tree_data`:

```python
def build_topology(
    sub_tree, start, tree_data, treefile, chrom_pos_dict, mut_data, purity=None
):
    """Recursively build one sample tree's topology dict.

    :param sub_tree:       clone id for this node
    :param start:          True only for the outermost call, forcing id 0
    :param tree_data:      one entry of summ.json's ``trees`` dict
    :param treefile:       per-tree JSON carrying ``mut_assignments``
    :param chrom_pos_dict: mutation-name to record mapping
    :param mut_data:       parsed mutation JSON carrying ``ssms``
    :param purity:         sample purity; unused until frequency normalisation
    :return: topology dict for this node and its descendants
    """
    if start:
        sub_tree = 0

    newsubtree = {
        "clone_id": int(sub_tree),
        "clone_mutations": [],
        "children": [],
        "X": 0,
        "x": 0,
        "new_x": 0,
    }

    if str(sub_tree) in tree_data["structure"]:
        for item in tree_data["structure"][str(sub_tree)]:

            child_dict = build_topology(
                item, False, tree_data, treefile, chrom_pos_dict, mut_data, purity
            )

            newsubtree["children"].append(child_dict)

        try:
            ssmli = []
            if start:
                pass
            else:
                for ssm in treefile["mut_assignments"][str(sub_tree)]["ssms"]:
                    ssmli.append(chrom_pos_dict[mut_data["ssms"][ssm]["name"]]["id"])
            newsubtree["clone_mutations"] = ssmli
            newsubtree["X"] = tree_data["populations"][str(sub_tree)][
                "cellular_prevalence"
            ][0]
            newsubtree["x"] = tree_data["populations"][str(sub_tree)][
                "cellular_prevalence"
            ][0]
            newsubtree["new_x"] = 0.0
        except Exception as e:
            print("Error in adding new subtree. Error not in base case**")
            print(sub_tree)
            print(e)
            pass

        return newsubtree

    else:
        # Base Case
        # make childrendict and return it
        ssmli = []

        for ssm in treefile["mut_assignments"][str(sub_tree)]["ssms"]:
            try:
                ssmli.append(chrom_pos_dict[mut_data["ssms"][ssm]["name"]]["id"])
            except Exception as e:
                print(
                    "Error in appending to mutation list. Error in base case appending ssm to ssmli"
                )
                print(e)
                # print(str(subTree))
                pass

        try:
            newsubtree["clone_mutations"] = ssmli
            newsubtree["X"] = tree_data["populations"][str(sub_tree)][
                "cellular_prevalence"
            ][0]
            newsubtree["x"] = tree_data["populations"][str(sub_tree)][
                "cellular_prevalence"
            ][0]
            newsubtree["new_x"] = 0.0
        except Exception as e:
            print("Error in adding new subtree. Error in base case")
            print(e)
            pass
        return newsubtree
```

- [ ] **Step 4: Update the call site**

In `main`, the `for tree in trees:` loop (originally line 292) becomes:

```python
    for tree in trees:

        inner_sample_tree_dict = {"topology": [], "score": trees[tree]["llh"]}
        with open("./" + args.tree_directory + "/" + str(tree) + ".json", "r") as f:
            # Load the JSON data into a dictionary
            treefile = json.load(f)

        bigtree = build_topology(
            tree, True, trees[tree], treefile, chrom_pos_dict, mut_data
        )

        inner_sample_tree_dict["topology"] = bigtree

        outer_dict["sample_trees"].append(inner_sample_tree_dict)
```

- [ ] **Step 5: Run the full suite**

Run: `pytest -v`
Expected: PASS — 8 tests in `test_generate_input.py`, plus all pre-existing tests.

- [ ] **Step 6: Confirm no leftover references**

Run: `grep -n "makeChild" src/neoantigen_utils/generate_input.py`
Expected: no output.

- [ ] **Step 7: Commit**

```bash
git add src/neoantigen_utils/generate_input.py tests/test_generate_input.py
git commit -m "refactor: lift makeChild out of main as build_topology"
```

---

### Task 3: Normalized X and exclusive x

**Files:**
- Modify: `src/neoantigen_utils/generate_input.py` — the two frequency blocks inside `build_topology`, the `build_topology` call site in `main`, and `VERSION` at line 26
- Test: `tests/test_generate_input.py` (extend)

**Interfaces:**
- Consumes: `compute_purity` (Task 1), `build_topology` and the
  `make_builder_inputs` / `walk_nodes` test helpers (Task 2).
- Produces: `purity` becomes a required-in-practice argument to
  `build_topology`; the default stays `None` so Task 2's signature is unchanged,
  but `main` always passes it.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_generate_input.py`:

```python
def assert_frequency_invariants(topology):
    """Check the clone-frequency contract on a built topology.

    - root X is exactly 1.0
    - every node: x == X - sum(children X)
    - leaves: x == X
    """
    assert topology["clone_id"] == 0
    assert topology["X"] == 1.0

    for node in walk_nodes(topology):
        children = node.get("children", [])
        expected_x = node["X"] - sum(c["X"] for c in children)
        assert node["x"] == pytest.approx(expected_x), (
            f"clone {node['clone_id']}: x={node['x']} != X - sum(children X)="
            f"{expected_x}"
        )
        if not children:
            assert node["x"] == pytest.approx(node["X"])


class TestFrequencyNormalization:
    def build(self):
        tree_data, treefile, chrom_pos_dict, mut_data = make_builder_inputs()
        purity = compute_purity(tree_data)
        return build_topology(
            "0", True, tree_data, treefile, chrom_pos_dict, mut_data, purity
        )

    def test_invariants_hold(self):
        assert_frequency_invariants(self.build())

    def test_expected_values(self):
        by_id = {n["clone_id"]: n for n in walk_nodes(self.build())}

        # purity == 0.8 + 0.2 == 1.0, so X equals the raw prevalence here
        assert by_id[0]["X"] == 1.0
        assert by_id[1]["X"] == pytest.approx(0.8)
        assert by_id[2]["X"] == pytest.approx(0.5)
        assert by_id[3]["X"] == pytest.approx(0.2)

        assert by_id[0]["x"] == pytest.approx(0.0)  # 1.0 - (0.8 + 0.2)
        assert by_id[1]["x"] == pytest.approx(0.3)  # 0.8 - 0.5
        assert by_id[2]["x"] == pytest.approx(0.5)  # leaf
        assert by_id[3]["x"] == pytest.approx(0.2)  # leaf

    def test_root_children_X_sum_to_one(self):
        topology = self.build()
        assert sum(c["X"] for c in topology["children"]) == pytest.approx(1.0)

    def test_normalization_divides_by_purity(self):
        # Halve every tumour prevalence: purity halves too, so normalised X is
        # unchanged. This is what distinguishes normalised X from raw prevalence.
        tree_data = make_tree_data(
            FIXTURE_STRUCTURE,
            {"0": 1.0, "1": 0.4, "2": 0.25, "3": 0.1},
        )
        _, treefile, chrom_pos_dict, mut_data = make_builder_inputs()
        purity = compute_purity(tree_data)
        assert purity == pytest.approx(0.5)

        topology = build_topology(
            "0", True, tree_data, treefile, chrom_pos_dict, mut_data, purity
        )
        by_id = {n["clone_id"]: n for n in walk_nodes(topology)}

        assert by_id[1]["X"] == pytest.approx(0.8)
        assert by_id[2]["X"] == pytest.approx(0.5)
        assert by_id[3]["X"] == pytest.approx(0.2)
        assert_frequency_invariants(topology)
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_generate_input.py::TestFrequencyNormalization -v`
Expected: FAIL — root `X` is the raw `1.0` by luck, but `test_expected_values`
fails on clone 1 (`x` is `0.8`, not `0.3`) and `test_normalization_divides_by_purity`
fails on clone 1 (`X` is `0.4`, not `0.8`).

- [ ] **Step 3: Replace the frequency assignment in the recursive branch**

Inside `build_topology`, in the first `try` block, replace:

```python
            newsubtree["X"] = tree_data["populations"][str(sub_tree)][
                "cellular_prevalence"
            ][0]
            newsubtree["x"] = tree_data["populations"][str(sub_tree)][
                "cellular_prevalence"
            ][0]
            newsubtree["new_x"] = 0.0
```

with:

```python
            newsubtree["X"] = (
                1.0
                if int(sub_tree) == 0
                else tree_data["populations"][str(sub_tree)]["cellular_prevalence"][0]
                / purity
            )
            newsubtree["x"] = newsubtree["X"] - sum(
                c["X"] for c in newsubtree["children"]
            )
            newsubtree["new_x"] = 0.0
```

Leave the surrounding `try`/`except` exactly as it is.

- [ ] **Step 4: Replace the frequency assignment in the base case**

Apply the identical replacement in the second `try` block (the base case). In
the base case `newsubtree["children"]` is empty, so the `sum(...)` is `0` and
`x == X`, which is the correct leaf behaviour:

```python
            newsubtree["X"] = (
                1.0
                if int(sub_tree) == 0
                else tree_data["populations"][str(sub_tree)]["cellular_prevalence"][0]
                / purity
            )
            newsubtree["x"] = newsubtree["X"] - sum(
                c["X"] for c in newsubtree["children"]
            )
            newsubtree["new_x"] = 0.0
```

- [ ] **Step 5: Pass purity from the call site**

In `main`'s `for tree in trees:` loop, compute purity before building and pass
it in:

```python
        try:
            purity = compute_purity(trees[tree])
        except ValueError as e:
            raise ValueError(f"Tree {tree}: {e}") from e

        bigtree = build_topology(
            tree, True, trees[tree], treefile, chrom_pos_dict, mut_data, purity
        )
```

- [ ] **Step 6: Bump the version marker**

`src/neoantigen_utils/generate_input.py:26` currently reads `VERSION = 1.9`.
Change it to:

```python
VERSION = 2.0
```

This is an output-format change, so it warrants the bump.

- [ ] **Step 7: Run the tests**

Run: `pytest tests/test_generate_input.py -v`
Expected: PASS, 12 tests.

Note: Task 2's `TestBuildTopologyStructure` tests pass `purity=1.0` and assert
only on shape, mutations, `new_x`, and the absence of `tilde_x`, so they remain
valid unchanged.

- [ ] **Step 8: Verify the test catches a regression**

Temporarily change the `x` line in the base case to
`newsubtree["x"] = newsubtree["X"]`.

Run: `pytest tests/test_generate_input.py::TestFrequencyNormalization -v`
Expected: still PASS — the base case only ever has empty children, so this edit
is a no-op there. Now make the same change in the **recursive branch** instead.
Expected: FAIL on `test_expected_values` and `test_invariants_hold` (clone 1).

Restore the correct implementation and re-run. Expected: PASS.

- [ ] **Step 9: Run the full suite**

Run: `pytest -v`
Expected: all PASS.

- [ ] **Step 10: Commit**

```bash
git add src/neoantigen_utils/generate_input.py tests/test_generate_input.py
git commit -m "fix: normalize X by purity and compute exclusive x per clone"
```

---

### Task 4: Document the output contract

**Files:**
- Modify: `README.md`

**Interfaces:**
- Consumes: nothing. Documentation only.

- [ ] **Step 1: Add a clone-frequency section to the README**

Find the section documenting `generate_input.py` and add:

```markdown
#### Clone frequency fields

Each node of `sample_trees[].topology` carries two frequencies:

- **`X`** — inclusive clone frequency: the PhyloWGS cellular prevalence divided
  by sample purity. The root clone is pinned to `1.0`, so the root's children
  sum to exactly `1.0`.
- **`x`** — exclusive clone frequency: `X` minus the sum of the children's `X`.
  For a leaf, `x == X`. For the root, `x == 0.0`.

Purity is the sum of the raw cellular prevalences of the root clone's children.
A tree whose purity is `0` is an error, not a tree of zeros.

`new_x` is emitted as `0.0` and populated downstream.

`tilde_x`, present in some NeoantigenEditing inputs, is **not** emitted. It is
only consumed by the paired primary-versus-recurrent longitudinal analysis, and
the rule determining which clones it zeroes is not available. See
`docs/superpowers/specs/2026-07-30-clone-frequency-normalization-design.md`.

The sample index is fixed at `cellular_prevalence[0]`; multi-sample
reconstructions are not supported.
```

- [ ] **Step 2: Commit**

```bash
git add README.md
git commit -m "docs: document clone frequency output contract"
```

---

## Verification

After all tasks:

- [ ] `pytest -v` passes from the repo root
- [ ] `git diff main --stat` touches only `generate_input.py`, `tests/test_generate_input.py`, `README.md`, and the spec/plan docs
- [ ] `grep -rn 'tilde_x' src/` returns nothing
- [ ] `grep -n 'makeChild' src/` returns nothing
- [ ] `grep -n 'VERSION' src/neoantigen_utils/generate_input.py` shows `2.0`

## Manual Sanity Check (optional, requires a real PhyloWGS run)

On real `summ.json` output, the load-bearing property is that each tree's root
children `X` sum to `1.0` and every node satisfies `x == X - sum(children X)` —
the same two checks that hold across all 57 nodes of the approved
NeoantigenEditing reference file.
