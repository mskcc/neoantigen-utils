"""Smoke tests for the accessory scripts folded into the package.

These verify the modules import cleanly and their console entry points resolve
(catching packaging/entry-point regressions). Fuller functional tests for
convertannotjson / format_netmhcpan are follow-up work.
"""

import subprocess


def test_convertannotjson_imports_and_has_entry():
    from neoantigen_utils import convertannotjson as m

    assert callable(m.main)
    assert callable(m.process_json_file)


def test_generate_input_imports_and_has_console_entry():
    # Regression guard: the module must import without pyensembl/mutalyzer (heavy
    # deps are tolerant), and expose the zero-arg ``console`` the entry point uses.
    from neoantigen_utils import generate_input as gi

    assert callable(gi.main)
    assert callable(gi.console)


def test_generate_mut_fasta_imports_and_has_entry():
    from neoantigen_utils import generate_mut_fasta as gmf

    assert callable(gmf.main)


def test_format_netmhcpan_imports_and_has_console_entry():
    from neoantigen_utils import format_netmhcpan_output as m

    assert callable(m.console)
    assert callable(m.parse_args)
    assert callable(m.netMHCpan_out_reformat)


def test_hla_string_imports_and_has_console_entry():
    # Fuller functional tests live in test_hla_string.py.
    from neoantigen_utils import hla_string as hs

    assert callable(hs.main)
    assert callable(hs.console)
    assert callable(hs.generate_hla_string)
