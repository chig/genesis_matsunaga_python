"""File path validation utilities for GENESIS Python interface."""

import glob
import os
from pathlib import Path
from typing import Dict, List, Optional, Union

from .exceptions import GenesisValidationError


def validate_file_exists(
    path: Optional[Union[str, Path, List[Union[str, Path]]]],
    name: str,
    required: bool = True
) -> None:
    """Validate that a file or list of files exists.

    Args:
        path: File path or list of file paths to validate
        name: Parameter name for error messages
        required: If True, raise error when path is None

    Raises:
        GenesisValidationError: If file does not exist or path is None when required
    """
    if path is None:
        if required:
            raise GenesisValidationError(
                f"{name} is required but was not provided"
            )
        return

    # Handle list of files (e.g., parfile, strfile)
    if isinstance(path, list):
        for i, p in enumerate(path):
            validate_file_exists(p, f"{name}[{i}]", required=True)
        return

    path = Path(path)
    if not path.exists():
        raise GenesisValidationError(
            f"{name}: File not found: {path}\n"
            f"Please check the file path and ensure the file exists."
        )
    if not path.is_file():
        raise GenesisValidationError(
            f"{name}: Path exists but is not a file: {path}"
        )


def validate_files_exist(
    files: Dict[str, Optional[Union[str, Path]]],
    required: Optional[List[str]] = None
) -> None:
    """Validate multiple files at once.

    Args:
        files: Dict mapping parameter names to file paths
        required: List of parameter names that are required (default: none required)
    """
    required = required or []
    for name, path in files.items():
        is_required = name in required
        validate_file_exists(path, name, required=is_required)


def validate_file_readable(path: Union[str, Path], name: str) -> None:
    """Validate that a file exists and is readable.

    Args:
        path: File path to validate
        name: Parameter name for error messages

    Raises:
        GenesisValidationError: If file is not readable
    """
    validate_file_exists(path, name, required=True)
    path = Path(path)
    if not os.access(path, os.R_OK):
        raise GenesisValidationError(
            f"{name}: File is not readable: {path}"
        )


def validate_file_pattern(
    pattern: Optional[str],
    name: str,
    required: bool = True
) -> None:
    """Validate a file pattern (e.g., 'file_{}.dat').

    Checks that at least one file matching the pattern exists.

    Args:
        pattern: File pattern with {} placeholder
        name: Parameter name for error messages
        required: If True, raise error when pattern is None

    Raises:
        GenesisValidationError: If no matching files found

    Note:
        This function checks if files exist at the pattern location.
        However, for patterns with format specifiers (like {}, {:d}, {0}),
        the actual file existence is checked by trying index 1 as a test.
        If no files exist yet (e.g., they will be created by the analysis),
        this validation is skipped.
    """
    if pattern is None:
        if required:
            raise GenesisValidationError(
                f"{name} is required but was not provided"
            )
        return

    import re

    # Check if pattern contains format placeholder
    # Matches: {}, {0}, {1}, {:d}, {:02d}, {0:d}, etc.
    format_pattern = r'\{[^}]*\}'
    if re.search(format_pattern, pattern):
        # Convert format placeholders to glob wildcards
        base_pattern = re.sub(format_pattern, '*', pattern)
        matches = glob.glob(base_pattern)
        if not matches:
            # No files found with glob - try with index 1 as fallback
            # This handles cases like {:d} that need actual formatting
            try:
                test_path = pattern.format(1)
                if not os.path.exists(test_path):
                    # Files don't exist yet - skip validation
                    # (they may be created by the analysis function)
                    return
            except (ValueError, KeyError, IndexError):
                # Pattern cannot be formatted with integer - skip validation
                return
    else:
        # Not a pattern, just a regular file
        validate_file_exists(pattern, name, required=True)


def validate_topology_combination(
    psffile: Optional[str] = None,
    prmtopfile: Optional[str] = None,
    grotopfile: Optional[str] = None,
) -> str:
    """Validate that at least one topology format is provided.

    Args:
        psffile: CHARMM PSF file path
        prmtopfile: AMBER PRMTOP file path
        grotopfile: GROMACS topology file path

    Returns:
        The detected format: 'CHARMM', 'AMBER', or 'GROMACS'

    Raises:
        GenesisValidationError: If no topology file is provided
    """
    has_charmm = psffile is not None
    has_amber = prmtopfile is not None
    has_gromacs = grotopfile is not None

    if not (has_charmm or has_amber or has_gromacs):
        raise GenesisValidationError(
            "At least one topology file is required: "
            "psffile (CHARMM), prmtopfile (AMBER), or grotopfile (GROMACS)"
        )

    # Return the detected format
    if has_charmm:
        return 'CHARMM'
    elif has_amber:
        return 'AMBER'
    else:
        return 'GROMACS'


def validate_trajectory_files(
    dcdfile: Optional[str] = None,
    trjfile: Optional[str] = None,
    name: str = "trajectory"
) -> None:
    """Validate trajectory file(s) exist.

    Args:
        dcdfile: DCD trajectory file path
        trjfile: Generic trajectory file path
        name: Parameter name for error messages

    Raises:
        GenesisValidationError: If no valid trajectory file found
    """
    if dcdfile is not None:
        validate_file_exists(dcdfile, f"{name} (dcdfile)")
    elif trjfile is not None:
        validate_file_exists(trjfile, f"{name} (trjfile)")
    # If neither is provided, that's okay - validation should happen
    # at the function level where we know if trajectory is required
