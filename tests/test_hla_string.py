"""Tests for neoantigen_utils.hla_string.

Fixture text mirrors the real POLYSOLVER / HLA-HD samples used to regression-test
the (now-retired) bash generateHLAString.sh in the modules repo, including the
three-digit-second-field case that was silently mis-scored before the fix.
"""

import sys

import pytest

from neoantigen_utils import hla_string

# Real POLYSOLVER winners.hla.txt shape: "HLA-<gene>\t<allele1>\t<allele2>".
POLYSOLVER_TWO_DIGIT = (
    "HLA-A\thla_a_24_02_01_01\thla_a_24_02_01_01\n"
    "HLA-B\thla_b_39_01_01_02\thla_b_39_01_01_02\n"
    "HLA-C\thla_c_07_01_05\thla_c_06_02_01_01\n"
)

# A three-digit second field (hla_b_18_177, hla_c_04_320, hla_c_07_348) used to be
# truncated by `cut -c 1-11` to HLA-B18:17 / HLA-C04:32 / HLA-C07:34 -- silently
# scoring the sample against an allele the patient may not carry.
POLYSOLVER_THREE_DIGIT = (
    "HLA-A\thla_a_02_01_01\thla_a_02_01_01\n"
    "HLA-B\thla_b_08_01_01\thla_b_18_177\n"
    "HLA-C\thla_c_04_320\thla_c_07_348\n"
)

# HLA-HD *_final.result.txt: bare locus label, '*'/':'-named alleles, "Not typed"
# for untyped loci, and possibly a second equally-scoring pair per line.
HLAHD_RESULT = (
    "A\tHLA-A*02:01:01\tHLA-A*24:02:01\n"
    "B\tHLA-B*07:02:01\tHLA-B*08:01:01\tHLA-B*07:05:01\tHLA-B*08:01:01\n"
    "C\tHLA-C*07:01:01\tHLA-C*07:02:01\n"
    "DRB1\tNot typed\tNot typed\n"
    "DQB1\tNot typed\tNot typed\n"
)


def test_polysolver_two_digit_unchanged():
    assert hla_string.generate_hla_string(POLYSOLVER_TWO_DIGIT) == (
        "HLA-A24:02,HLA-A24:02,HLA-B39:01,HLA-B39:01,HLA-C07:01,HLA-C06:02"
    )


def test_polysolver_three_digit_second_field_not_truncated():
    # The whole point of the fix: a three-digit second field must survive intact.
    assert hla_string.generate_hla_string(POLYSOLVER_THREE_DIGIT) == (
        "HLA-A02:01,HLA-A02:01,HLA-B08:01,HLA-B18:177,HLA-C04:320,HLA-C07:348"
    )


def test_polysolver_duplicate_alleles_preserved_for_homozygous_samples():
    result = hla_string.generate_hla_string(POLYSOLVER_TWO_DIGIT)
    assert result.split(",").count("HLA-A24:02") == 2


def test_polysolver_unparseable_entry_is_skipped_and_warned(capsys):
    text = "HLA-A\thla_a_02_01_01\tgarbage\n"
    result = hla_string.generate_hla_string(text)
    assert result == "HLA-A02:01"
    assert "garbage" in capsys.readouterr().err


def test_hlahd_output_shape_matches_polysolver():
    assert hla_string.generate_hla_string(HLAHD_RESULT) == (
        "HLA-A02:01,HLA-A24:02,HLA-B07:02,HLA-B08:01,HLA-C07:01,HLA-C07:02"
    )


def test_hlahd_class_ii_loci_are_dropped():
    result = hla_string.generate_hla_string(HLAHD_RESULT)
    assert "DRB" not in result and "DQB" not in result


def test_hlahd_only_first_pair_per_locus_is_emitted():
    # The B line lists a second, equally-scoring pair; only the first is used.
    result = hla_string.generate_hla_string(HLAHD_RESULT)
    assert sum(a.startswith("HLA-B") for a in result.split(",")) == 2


def test_hlahd_not_typed_is_dropped_silently(capsys):
    result = hla_string.generate_hla_string(HLAHD_RESULT)
    assert "Not typed" not in result
    assert "Not typed" not in capsys.readouterr().err


def test_hlahd_class_ii_locus_is_dropped_silently(capsys):
    # DRB1/DQB1/etc are an expected, routine part of real HLA-HD output --
    # not a parse failure -- so dropping them must stay silent.
    hla_string.generate_hla_string(HLAHD_RESULT)
    assert capsys.readouterr().err == ""


def test_hlahd_malformed_line_is_warned_not_silently_dropped(capsys):
    # A line with no allele columns at all can't be a real HLA-HD row (every
    # real row -- even a fully "Not typed" one -- has at least one allele
    # field). Silently treating it the same as an expected class II skip
    # would hide a misdetected format or a corrupted file with no diagnostic
    # trail.
    text = "A\tHLA-A*02:01:01\tHLA-A*24:02:01\nXYZ\n"
    result = hla_string.generate_hla_string(text)
    assert result == "HLA-A02:01,HLA-A24:02"
    assert "XYZ" in capsys.readouterr().err


def test_empty_input_raises():
    with pytest.raises(ValueError):
        hla_string.generate_hla_string("")


def test_label_only_input_raises():
    with pytest.raises(ValueError):
        hla_string.generate_hla_string("HLA-A\nHLA-B\nHLA-C\n")


def test_main_writes_result_to_stdout_and_exits_zero(tmp_path, capsys):
    input_file = tmp_path / "winners.hla.txt"
    input_file.write_text(POLYSOLVER_TWO_DIGIT)

    hla_string.main(["-f", str(input_file)])

    assert capsys.readouterr().out.strip() == (
        "HLA-A24:02,HLA-A24:02,HLA-B39:01,HLA-B39:01,HLA-C07:01,HLA-C06:02"
    )


def test_main_exits_one_on_unparseable_file(tmp_path):
    input_file = tmp_path / "empty.txt"
    input_file.write_text("")

    with pytest.raises(SystemExit) as exc:
        hla_string.main(["-f", str(input_file)])
    assert exc.value.code == 1


def test_main_dash_f_empty_string_is_an_error_not_usage(capsys):
    # -f "" is an explicitly-passed, unusable path -- not the same as omitting
    # -f entirely. Treating them alike would let a pipeline that hands in an
    # unresolved/empty variable "succeed" with no alleles, or feed usage text
    # into netMHCpan's -a argument.
    with pytest.raises(SystemExit) as exc:
        hla_string.main(["-f", ""])
    assert exc.value.code == 1
    assert "USAGE" not in capsys.readouterr().out


def test_main_missing_file_exits_one_with_structured_error(capsys):
    with pytest.raises(SystemExit) as exc:
        hla_string.main(["-f", "/no/such/file.txt"])
    assert exc.value.code == 1
    assert "ERROR" in capsys.readouterr().err


def test_console_prints_version_and_exits_zero(monkeypatch, capsys):
    monkeypatch.setattr(sys, "argv", ["generateHLAString.sh", "-v"])
    with pytest.raises(SystemExit) as exc:
        hla_string.console()
    assert exc.value.code == 0
    assert capsys.readouterr().out.strip() == hla_string.VERSION
