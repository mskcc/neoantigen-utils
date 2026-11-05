"""find_most_similar_string must return exactly what the original (full
alignment, score taken from alignments[0][2]) implementation returned."""
import random

import pytest

pairwise2 = pytest.importorskip("Bio.pairwise2")

from neoantigen_utils import generate_input as gi

AA = "ACDEFGHIKLMNPQRSTVWY"


def _reference(target, strings):
    # Verbatim pre-change logic, scoring via the full alignment list.
    max_score, max_score2 = -1, -2
    best = best2 = first_same = None
    first_same_score = -1
    for s in strings:
        if len(s) == len(target):
            score = pairwise2.align.globalxx(target, s)[0][2]
            if score > max_score2:
                if score > max_score:
                    max_score2, best2 = max_score, best
                    max_score, best = score, s
                else:
                    max_score2, best2 = score, s
            if target[0] == s[0] and score > first_same_score:
                first_same, first_same_score = s, score
    return best, best2, first_same, first_same_score, max_score


def _mutate(p, rng, n):
    p = list(p)
    for _ in range(n):
        p[rng.randrange(len(p))] = rng.choice(AA)
    return "".join(p)


def _cases():
    rng = random.Random(0)
    for _ in range(150):
        k = rng.choice([8, 9, 10, 11])
        target = "".join(rng.choice(AA) for _ in range(k))
        # WT-like windows: near-copies (many tied scores), randoms, other lengths
        cands = [_mutate(target, rng, rng.randint(1, 4)) for _ in range(6)]
        cands += ["".join(rng.choice(AA) for _ in range(k)) for _ in range(3)]
        cands += ["".join(rng.choice(AA) for _ in range(k + 1))]
        cands += [cands[0]]  # duplicate -> exact tie
        rng.shuffle(cands)
        yield target, cands
    yield "ACDEFGHIK", []
    yield "ACDEFGHIK", ["ACDEFGHIKL"]  # no same-length candidate


@pytest.mark.parametrize("target,cands", list(_cases()))
def test_matches_reference(target, cands):
    got = gi.find_most_similar_string(target, cands)
    assert got == _reference(target, cands)
    assert all(type(a) is type(b) for a, b in zip(got, _reference(target, cands)))
