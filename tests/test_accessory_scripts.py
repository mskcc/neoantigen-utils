"""Smoke tests for the accessory scripts folded into the package.

These verify the modules import cleanly and their console entry points resolve
(catching packaging/entry-point regressions). Fuller functional tests for
convertannotjson / format_netmhcpan are follow-up work.
"""

import subprocess

import pytest


def test_convertannotjson_imports_and_has_entry():
    from neoantigen_utils import convertannotjson as m

    assert callable(m.main)
    assert callable(m.process_json_file)


def test_format_netmhcpan_imports_and_has_console_entry():
    from neoantigen_utils import format_netmhcpan_output as m

    assert callable(m.console)
    assert callable(m.parse_args)
    assert callable(m.netMHCpan_out_reformat)


def test_hla_string_shim_execs_bundled_script_version():
    # The shim forwards to the bundled generateHLAString.sh; ``-v`` prints the
    # version and exits 0.
    from neoantigen_utils import hla_string

    with pytest.raises(SystemExit) as exc:
        import sys

        old = sys.argv
        sys.argv = ["generateHLAString.sh", "-v"]
        try:
            hla_string.console()
        finally:
            sys.argv = old
    assert exc.value.code == 0


def test_bundled_hla_script_is_present():
    from importlib import resources

    path = resources.files("neoantigen_utils.data") / "generateHLAString.sh"
    assert path.is_file()
