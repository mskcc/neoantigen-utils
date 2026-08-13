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
    """Sniff whether ``text`` is HLA-HD (``HLA-<gene>*...``) rather than
    POLYSOLVER output, by looking for HLA-HD's asterisk-named alleles
    anywhere in the file."""
    return bool(re.search(r"HLA-[A-Za-z0-9]*\*", text))


def parse_polysolver_allele(item: str):
    """Parse one raw POLYSOLVER allele field, e.g. ``hla_a_02_01_01``.

    Shared by :func:`parse_polysolver` and
    :func:`neoantigen_utils.generate_input.convert_polysolver_hla`, which
    format the parsed fields differently for their own consumers -- keeping
    the parsing itself in one place means a future POLYSOLVER format fix
    (e.g. the 3-digit-field truncation bug described in this module's
    docstring) only needs to be made once.

    :param item: one raw allele field, not the ``HLA-<gene>`` label token
    :return: ``(gene, field1, field2)``, all upper-cased, or ``None`` (after
             printing a warning to stderr) if ``item`` doesn't match the
             expected ``hla_<gene>_<field1>_<field2>[_...]`` shape
    """
    parts = item.upper().split("_")
    if len(parts) < 4 or parts[0] != "HLA" or not parts[1] or not parts[2] or not parts[3]:
        print(f"WARN: skipping unparseable HLA entry '{item}'", file=sys.stderr)
        return None
    return parts[1], parts[2], parts[3]


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

        parsed = parse_polysolver_allele(item)
        if parsed is None:
            continue

        gene, field1, field2 = parsed
        alleles.append(f"HLA-{gene}{field1}:{field2}")

    return alleles


def parse_hlahd(text: str) -> list[str]:
    """Parse HLA-HD *_final.result.txt content into two-field allele names."""
    alleles = []
    lines_seen = 0
    class_i_lines_seen = 0
    for line in text.replace("\r", "").split("\n"):
        if not line:
            continue
        lines_seen += 1

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
            # Expected: HLA-HD reports every locus (DRB1, DQB1, ...), and only
            # A/B/C are class I. Not warned per-line -- that would fire on
            # every well-formed file -- see the aggregate check below instead.
            continue
        class_i_lines_seen += 1

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

    if lines_seen and not class_i_lines_seen:
        # Every line was dropped by the locus filter -- most likely this file
        # was misdetected as HLA-HD (or genuinely has no A/B/C rows), and the
        # per-line skip above would otherwise leave no trace of why.
        print(
            f"WARN: no HLA-A/B/C locus lines found among {lines_seen} line(s) "
            "of HLA-HD input",
            file=sys.stderr,
        )

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
    """Parse CLI arguments.

    :param argv: argument list, or ``None`` to use ``sys.argv``
    :return: parsed ``Namespace``, with ``file`` set to ``None`` if ``-f``
             was omitted entirely (an explicitly empty ``-f ""`` is left as
             ``""``, distinct from omission, so :func:`main` can fail loudly
             on it instead of silently printing usage)
    """
    parser = argparse.ArgumentParser(
        prog="generateHLAString.sh",
        description="Build the netMHCpan -a allele string from an HLA typing result.",
    )
    parser.add_argument("-f", dest="file", default=None, help="HLA typing result file")
    parser.add_argument("-v", "--version", action="version", version=VERSION)
    return parser.parse_args(argv)


def main(argv=None):
    """CLI entry point: parse args, resolve the HLA string, print it or exit(1).

    :param argv: argument list, or ``None`` to use ``sys.argv``
    """
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
    """Zero-argument console entry point (parses ``sys.argv`` then runs ``main``)."""
    main()


if __name__ == "__main__":
    console()
