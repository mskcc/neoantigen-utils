"""Unit tests for the NMDetective-B (Lindeboom et al. 2019) implementation in
``generate_input.py``.

The module under test imports several heavy third-party libraries at import
time (Bio, pyensembl, ...). Mirroring the import-stubbing pattern used
elsewhere in this repo, we inject lightweight stubs into ``sys.modules`` for any
dependency that may be absent (notably ``pyensembl``) before loading the module
by file path. The NMD logic we test does not use those libraries directly.
"""

import csv
import importlib.util
import sys
import types
from pathlib import Path

import pytest


# --------------------------------------------------------------------------- #
# Module loading with dependency stubbing
# --------------------------------------------------------------------------- #
def _install_stub(name):
    """Register a bare module (and parent packages) in sys.modules if missing."""
    if name in sys.modules:
        return
    parts = name.split(".")
    for i in range(1, len(parts) + 1):
        sub = ".".join(parts[:i])
        if sub not in sys.modules:
            sys.modules[sub] = types.ModuleType(sub)


def _load_generate_input():
    # pyensembl is not required for the pure NMD logic; stub it unconditionally
    # so the module imports on hosts without pyensembl / a GRCh37 cache.
    for mod in ("pyensembl", "pyensembl.genome"):
        _install_stub(mod)
    if not hasattr(sys.modules["pyensembl"], "EnsemblRelease"):
        sys.modules["pyensembl"].EnsemblRelease = object
    if not hasattr(sys.modules["pyensembl.genome"], "Genome"):
        sys.modules["pyensembl.genome"].Genome = object

    # Bio / numpy are usually present, but stub Bio pieces defensively.
    try:
        import Bio  # noqa: F401
        import Bio.pairwise2  # noqa: F401
    except Exception:
        _install_stub("Bio")
        _install_stub("Bio.pairwise2")
        sys.modules["Bio.pairwise2"].pairwise2 = types.SimpleNamespace()
        sys.modules["Bio.pairwise2"].format_alignment = lambda *a, **k: ""
        sys.modules["Bio"].pairwise2 = sys.modules["Bio.pairwise2"]

    here = Path(__file__).resolve().parent
    script = (
        here.parent / "src" / "neoantigen_utils" / "generate_input.py"
    )
    spec = importlib.util.spec_from_file_location("generate_input", script)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


gi = _load_generate_input()


FIXTURE = Path(__file__).resolve().parent / "nmd_table2_fixture.csv"

# Map the fixture's label vocabulary onto the module's exact return strings.
LABEL_MAP = {
    "Last exon": "Last Exon",
    "Start-proximal": "Start-proximal",
    "Long exon": "Long Exon",
    "50nt Rule": "50nt Rule",
    "Trigger NMD": "Trigger NMD",
}


def _load_fixture_rows():
    with open(FIXTURE, newline="") as fh:
        return list(csv.DictReader(fh))


# --------------------------------------------------------------------------- #
# 1. Ground-truth table (paper Suppl. Table 2) drives the pure classifier
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("row", _load_fixture_rows(), ids=lambda r: f"{r['gene']}_{r['pos']}")
def test_classify_nmd_matches_paper_table2(row):
    is_last_exon = row["lastExon"].strip().lower() == "true"
    dist_from_start = int(row["distToStart"])
    exon_len = int(row["exonLength"])
    within_50 = row["rule50bp"].strip().lower() == "true"

    expected = LABEL_MAP[row["expected_nmd"].strip()]
    result = gi.classify_nmd(is_last_exon, dist_from_start, exon_len, within_50)
    assert result == expected


def test_fixture_covers_both_strands_and_two_labels():
    rows = _load_fixture_rows()
    strands = {r["strand"] for r in rows}
    labels = {r["expected_nmd"].strip() for r in rows}
    assert strands == {"+", "-"}
    assert labels == {"Last exon", "Start-proximal"}
    assert len(rows) == 23


# --------------------------------------------------------------------------- #
# 2a. Pure decision-tree branch coverage (branches the table lacks)
# --------------------------------------------------------------------------- #
def test_classify_last_exon_wins_over_everything():
    # Even with a near-start, short, penultimate PTC, last-exon takes priority.
    assert gi.classify_nmd(True, 10, 100, True) == "Last Exon"


def test_classify_start_proximal_boundary():
    assert gi.classify_nmd(False, 149, 1000, False) == "Start-proximal"
    # exactly 150 is NOT start-proximal
    assert gi.classify_nmd(False, 150, 1000, False) == "Long Exon"


def test_classify_long_exon():
    assert gi.classify_nmd(False, 300, 408, False) == "Long Exon"
    # 407 is not "long" (must be > 407)
    assert gi.classify_nmd(False, 300, 407, False) == "Trigger NMD"


def test_classify_50nt_rule():
    assert gi.classify_nmd(False, 300, 200, True) == "50nt Rule"


def test_classify_trigger_nmd():
    assert gi.classify_nmd(False, 300, 200, False) == "Trigger NMD"


# --------------------------------------------------------------------------- #
# 2b. Feature-extraction / end-to-end coverage on hand-built transcripts.
#     Each test targets a specific confirmed bug and would FAIL under the old
#     genomic-coordinate logic.
# --------------------------------------------------------------------------- #
def _classify_from_exons(exons, strand, mut_pos, mut_to_stop_dist, start_codon_offset=0):
    feats = gi.extract_nmd_features(
        exons, strand, mut_pos, mut_to_stop_dist, start_codon_offset
    )
    return feats, gi.classify_nmd(**feats)


def test_bug6_exon_length_off_by_one():
    # Inclusive genomic range 100..200 is 101 nt, not 100.
    assert gi._exon_len((100, 200)) == 101


def test_bug1_far_ptc_is_not_start_proximal_plus_strand():
    # Plus strand, 3 short exons of 100 nt each (offsets 0..99, 100..199, 200..299).
    # Mutation early in exon 1, PTC 220 nt downstream lands in exon 3 (the last
    # exon here) -> to isolate bug #1 we add a 4th exon so exon 3 is internal.
    exons = [(1000, 1099), (2000, 2099), (3000, 3099), (4000, 4099)]  # 4x100nt
    # mutation at spliced offset 10 (exon1). PTC 210 nt downstream => offset 220
    # => exon index 2 (third exon), internal, dist_from_start=220 (>150),
    # exon len 100 (<=407), not penultimate-within-50 => Trigger NMD.
    feats, label = _classify_from_exons(exons, "+", 1010, 210)
    assert feats["is_last_exon"] is False
    assert feats["dist_from_start_nt"] == 220
    assert label == "Trigger NMD"  # old sign-inverted logic returned Start-proximal


def test_bug5_minus_strand_spliced_offset():
    # Minus strand: transcript 5' end is the HIGH genomic coordinate.
    # exons ordered 5'->3': (2000,2099) then (1000,1099).
    exons = [(2000, 2099), (1000, 1099)]
    # genomic pos 2090 is 9 nt into the transcript (2099 -> offset0 ... 2090 -> offset9)
    assert gi._spliced_offset(exons, "-", 2090) == 9
    # genomic pos 1099 is start of exon2 => offset 100
    assert gi._spliced_offset(exons, "-", 1099) == 100
    # A plus-strand reading of the same coords would be wrong:
    assert gi._spliced_offset(exons, "+", 2090) == 90


def test_bug4_own_exon_ptc_anchored_to_mutation_not_exon_start():
    # Single large exon; PTC stays in the mutation's own exon and its distance
    # from the start codon must be measured from the mutation position + offset,
    # NOT from the exon start.
    exons = [(1000, 2000)]  # 1001 nt, last & only exon
    feats = gi.extract_nmd_features(exons, "+", 1500, 30, start_codon_offset=0)
    # mutation spliced offset = 500; PTC offset = 530; dist_from_start = 530.
    assert feats["dist_from_start_nt"] == 530
    assert feats["is_last_exon"] is True


def test_bug3_same_exon_by_index_not_length_equality():
    # Two exons of EQUAL length (100 nt). Mutation in exon 1; a short PTC must
    # stay in exon 1, decided by index, not by matching exon lengths.
    exons = [(1000, 1099), (5000, 5099)]  # both 100 nt
    feats = gi.extract_nmd_features(exons, "+", 1010, 20)  # PTC at offset 30
    # PTC is in exon index 0 (the last is index1) -> not last exon.
    assert feats["is_last_exon"] is False
    # exon length reported is exon 0's length (100), unambiguous.
    assert feats["ptc_exon_len_nt"] == 100


def test_bug2_50nt_rule_measures_distance_to_final_junction():
    # 3 exons; PTC lands in the penultimate exon within 50 nt of its 3' end
    # (the final exon-exon junction). Old code measured from the 5' start.
    # Exon layout (plus strand), spliced offsets:
    #   exon0: 0..99   (100 nt)
    #   exon1: 100..699 (600 nt)  <- penultimate; 3' end at offset 699
    #   exon2: 700..799 (100 nt)  <- last
    exons = [(1000, 1099), (2000, 2599), (3000, 3099)]
    # Put PTC at spliced offset 680 => 19 nt upstream of junction (offset 699).
    # mutation at offset 100 (start of exon1), PTC 580 downstream => 680.
    feats = gi.extract_nmd_features(exons, "+", 2000, 580)
    assert feats["is_last_exon"] is False
    assert feats["within_50nt_of_last_junction"] is True
    # But penultimate exon is 600 nt (>407): Long Exon wins per tree order.
    assert gi.classify_nmd(**feats) == "Long Exon"

    # Now a SHORT penultimate exon so the 50nt branch is actually reached.
    exons2 = [(1000, 1099), (2000, 2199), (3000, 3099)]  # penult = 200 nt
    # penult offsets 100..299, junction at 299. PTC at 280 => 19 nt upstream.
    feats2 = gi.extract_nmd_features(exons2, "+", 2000, 180)
    assert feats2["within_50nt_of_last_junction"] is True
    assert gi.classify_nmd(**feats2) == "50nt Rule"


def test_end_to_end_long_exon_branch():
    # PTC well past the start, in a long (>407) internal exon.
    exons = [(1000, 1099), (2000, 2999), (3000, 3099)]  # middle exon 1000 nt
    feats = gi.extract_nmd_features(exons, "+", 2010, 200)  # PTC offset ~ 110+200
    assert feats["ptc_exon_len_nt"] == 1000
    assert feats["is_last_exon"] is False
    assert gi.classify_nmd(**feats) == "Long Exon"


# --------------------------------------------------------------------------- #
# 3. Robustness: unlocalizable PTC -> "False", no exception.
# --------------------------------------------------------------------------- #
def test_extract_features_returns_none_for_intronic_mutation():
    exons = [(1000, 1099), (2000, 2099)]
    # genomic pos 1500 is in the intron between the two exons.
    assert gi.extract_nmd_features(exons, "+", 1500, 30) is None


def test_determine_nmd_returns_false_when_transcript_unresolvable():
    class _Ensembl:
        def transcript_by_id(self, tid):
            raise ValueError("no such transcript")

    result = gi.determine_NMD("1", 1500, 5, 0, _Ensembl(), transcriptID="ENST0000FAKE.1")
    assert result == "False"


def test_determine_nmd_returns_false_for_intronic_via_stub_ensembl():
    # Build a fake pyensembl-like transcript so determine_NMD runs end-to-end
    # and returns "False" for a mutation that is not inside any exon.
    Exon = lambda s, e: types.SimpleNamespace(start=s, end=e)

    class _T:
        strand = "+"
        exons = [Exon(1000, 1099), Exon(2000, 2099)]
        start_codon_spliced_offsets = [0]

    class _Ensembl:
        def transcript_by_id(self, tid):
            return _T()

    # pos 1500 is intronic -> unlocalizable -> "False", no exception.
    assert gi.determine_NMD("1", 1500, 5, 0, _Ensembl(), transcriptID="ENST1") == "False"


def test_determine_nmd_end_to_end_last_exon_via_stub_ensembl():
    Exon = lambda s, e: types.SimpleNamespace(start=s, end=e)

    class _T:
        strand = "+"
        # two exons; PTC in the last one
        exons = [Exon(1000, 1099), Exon(2000, 2999)]
        start_codon_spliced_offsets = [0]

    class _Ensembl:
        def transcript_by_id(self, tid):
            return _T()

    # mutation at 2010 (offset 110), short PTC stays in last exon.
    assert gi.determine_NMD("1", 2010, 5, 0, _Ensembl(), transcriptID="ENST1.3") == "Last Exon"


# --------------------------------------------------------------------------- #
# 4. Optional integration test: only runs if pyensembl + a GRCh37 cache exist.
# --------------------------------------------------------------------------- #
def test_optional_pyensembl_integration():
    pytest.importorskip("pyensembl", reason="pyensembl not installed")
    try:
        from pyensembl import EnsemblRelease
        ensembl = EnsemblRelease(75)  # GRCh37
        # Trigger a lookup; if the cache is missing this raises.
        ensembl.transcript_by_id("ENST00000218516")
    except Exception as exc:  # noqa: BLE001
        pytest.skip(f"GRCh37 pyensembl cache unavailable: {exc}")

    rows = _load_fixture_rows()
    checked = 0
    for row in rows[:3]:
        tid = row["transcript"].split(".")[0]
        try:
            t = ensembl.transcript_by_id(tid)
        except Exception:
            continue
        exons, strand, sco = gi._get_transcript_model(
            row["chr"].replace("chr", ""), int(row["pos"]), ensembl, tid
        )
        assert strand == row["strand"]
        checked += 1
    if checked == 0:
        pytest.skip("no fixture transcripts resolvable in this cache")
