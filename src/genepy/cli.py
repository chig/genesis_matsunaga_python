# SPDX-License-Identifier: LGPL-3.0-or-later
"""
CLI entry points for GENESIS commands.

This module provides Python entry points for bundled GENESIS binaries.
"""

import subprocess
import sys
from pathlib import Path


def get_bin_path(name: str) -> Path:
    """Get path to bundled binary.

    Args:
        name: Name of the binary (e.g., 'atdyn')

    Returns:
        Path to the binary executable
    """
    return Path(__file__).parent / "bin" / name


def run_atdyn():
    """Run atdyn command.

    Entry point for the atdyn molecular dynamics engine.
    All command line arguments are passed through to the binary.
    """
    binary = get_bin_path("atdyn")
    if not binary.exists():
        print(f"Error: atdyn binary not found at {binary}", file=sys.stderr)
        print("The binary may not be included in this installation.", file=sys.stderr)
        print("Please build GENESIS from source or install via conda.", file=sys.stderr)
        sys.exit(1)
    sys.exit(subprocess.call([str(binary)] + sys.argv[1:]))


# Note: spdyn requires MPI and is not included in PyPI package.
# Users who need spdyn should build GENESIS from source or use conda.
