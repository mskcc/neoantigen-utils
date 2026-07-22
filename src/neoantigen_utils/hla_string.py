"""Console shim for the bundled ``generateHLAString.sh``.

``pip`` installs Python console scripts, not shell scripts, so this exposes the
bundled shell script under its original command name (``generateHLAString.sh``)
by exec'ing it with bash and forwarding all arguments and the exit code.
"""

import subprocess
import sys
from importlib import resources


def console():
    with resources.as_file(
        resources.files("neoantigen_utils.data") / "generateHLAString.sh"
    ) as script:
        sys.exit(subprocess.call(["bash", str(script), *sys.argv[1:]]))


if __name__ == "__main__":
    console()
