"""Build the netMHCpan/netMHCstabpan `-a` allele argument from an HLA typing result.

Understands two input formats, auto-detected from file contents (no CLI flag; a
file comes from one caller and is never a mix of the two):

- POLYSOLVER ``winners.hla.txt``: one locus per line, a bare ``HLA-<gene>`` label
  followed by tab-separated allele fields, e.g. ``HLA-B\\thla_b_08_01_01\\thla_b_18_177``.
- HLA-HD ``*_final.result.txt``: one locus per line, a bare locus label (``A``,
  ``B``, ``C``, ...) followed by ``*``/``:``-named alleles, e.g.
  ``B\\tHLA-B*07:02:01\\tHLA-B*08:01:01``. Untyped loci carry the literal
  ``Not typed``, and a locus may list a second, equally-scoring pair.

This module is the single source of truth for this parsing: it used to be a bash
script (``cut -c 1-11``) that assumed the POLYSOLVER allele's second field was
always two digits. It is not -- ``hla_b_18_177`` truncated to ``HLA-B18:17``, a
real allele the patient may not carry, silently mis-scoring the sample.
"""

import argparse
import re
import sys

VERSION = "1.6.1"

_HLAHD_ALLELE_RE = re.compile(r"^HLA-([A-Za-z0-9]+)\*([0-9]+):([0-9]+)")
_HLAHD_CLASS_I_LOCI = {"A", "B", "C"}


def _is_hlahd_format(text: str) -> bool:
    """Sniff whether `text` is HLA-HD (vs. POLYSOLVER) by looking for `HLA-<gene>*`."""
    return bool(re.search(r"HLA-[A-Za-z0-9]*\*", text))


def parse_polysolver(text: str) -> list[str]:
    """Parse POLYSOLVER winners.hla.txt content into two-field allele names."""
    alleles = []
    # Flatten every line's tab-separated fields, then drop the bare "HLA-<gene>"
    # label tokens -- they use a hyphen, allele fields use an underscore
    # (hla_a_02_01_01), so a simple prefix check tells them apart.
    for item in text.replace("\r", "").replace("\t", "\n").split("\n"):
        if not item:
            continue
        if item.upper().startswith("HLA-"):
            continue

        parts = item.upper().split("_")
        if len(parts) < 4 or parts[0] != "HLA" or not parts[1] or not parts[2] or not parts[3]:
            print(f"WARN: skipping unparseable HLA entry '{item}'", file=sys.stderr)
            continue

        _prefix, gene, field1, field2 = parts[0], parts[1], parts[2], parts[3]
        alleles.append(f"HLA-{gene}{field1}:{field2}")

    return alleles


def parse_hlahd(text: str) -> list[str]:
    """Parse HLA-HD *_final.result.txt content into two-field allele names."""
    alleles = []
    for line in text.replace("\r", "").split("\n"):
        if not line:
            continue

        fields = line.split("\t")
        locus = fields[0]

        # Every real HLA-HD row has a locus label plus at least one allele
        # field -- even a fully untyped one ("A\tNot typed\tNot typed"). A row
        # with no allele column at all can't be a real row, so warn instead of
        # silently dropping it the same way an expected class II locus (DRB1,
        # DQB1, ...) is dropped.
        if len(fields) < 2:
            print(f"WARN: skipping unparseable HLA entry '{line}'", file=sys.stderr)
            continue

        if locus not in _HLAHD_CLASS_I_LOCI:
            continue

        # Only the first pair is taken. HLA-HD may report an alternative,
        # equally-scoring pair on the same line; emitting it would exceed
        # POLYSOLVER's two-per-locus cardinality and score the sample against
        # alleles it may not carry.
        for allele in fields[1:3]:
            if not allele or allele == "Not typed":
                continue

            match = _HLAHD_ALLELE_RE.match(allele)
            if not match:
                print(f"WARN: skipping unparseable HLA entry '{allele}'", file=sys.stderr)
                continue

            gene, field1, field2 = match.groups()
            alleles.append(f"HLA-{gene}{field1}:{field2}")

    return alleles


def generate_hla_string(text: str) -> str:
    """Return the comma-joined ``-a`` allele string for netMHCpan/netMHCstabpan.

    Raises ``ValueError`` if no alleles could be parsed from ``text``.
    """
    alleles = parse_hlahd(text) if _is_hlahd_format(text) else parse_polysolver(text)

    if not alleles:
        raise ValueError("no HLA alleles parsed from input")

    return ",".join(alleles)


def parse_args(argv=None):
    """Parse CLI args: `-f FILE` (required to do anything), `-v`/`--version`."""
    parser = argparse.ArgumentParser(
        prog="generateHLAString.sh",
        description="Build the netMHCpan -a allele string from an HLA typing result.",
    )
    parser.add_argument("-f", dest="file", required=False, help="HLA typing result file")
    parser.add_argument("-v", "--version", action="version", version=VERSION)
    return parser.parse_args(argv)


def main(argv=None):
    """Parse argv, print the netMHCpan `-a` allele string, and exit 1 on failure."""
    args = parse_args(argv)
    # `-f` omitted entirely (args.file is None) just prints usage. An explicitly
    # passed but empty path ("") is not the same thing -- it's an unusable
    # argument, not "no flag given" -- and must be treated as an error, not
    # silently accepted.
    if args.file is None:
        print("USAGE: generateHLASTRING.sh -f [HLA_FILE]")
        return

    try:
        with open(args.file) as f:
            text = f.read()
    except OSError as e:
        print(f"ERROR: could not read {args.file!r}: {e}", file=sys.stderr)
        sys.exit(1)

    try:
        result = generate_hla_string(text)
    except ValueError:
        print(f"ERROR: no HLA alleles parsed from {args.file}", file=sys.stderr)
        sys.exit(1)

    print(result)


def console():
    """Zero-argument console entry point (parses `sys.argv` then runs `main`)."""
    main()


if __name__ == "__main__":
    console()
