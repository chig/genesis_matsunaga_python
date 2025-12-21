"""
Validation utilities for GENESIS Python interface.
"""

import ctypes
from typing import Sequence, Union

import numpy as np

from .exceptions import (
    GenesisValidationError,
    GenesisMemoryError,
    GenesisOverflowError,
)


# Maximum safe size for array operations (2GB limit for 64-bit systems)
MAX_SAFE_SIZE = 2**31 - 1


def validate_size(size: Union[int, Sequence[int]], name: str = "size") -> int:
    """Validate array size and return total element count.

    Args:
        size: Single dimension or sequence of dimensions
        name: Name for error messages

    Returns:
        Total number of elements

    Raises:
        GenesisValidationError: If size is negative or empty
        GenesisOverflowError: If total size exceeds safe limits
    """
    if isinstance(size, (int, np.integer)):
        if size < 0:
            raise GenesisValidationError(
                f"{name} must be non-negative, got {size}"
            )
        if size > MAX_SAFE_SIZE:
            raise GenesisOverflowError(
                f"{name} ({size}) exceeds maximum safe size ({MAX_SAFE_SIZE})"
            )
        return int(size)

    # It's a sequence
    if len(size) == 0:
        raise GenesisValidationError(f"{name} cannot be empty sequence")

    total = 1
    for i, dim in enumerate(size):
        if dim < 0:
            raise GenesisValidationError(
                f"{name}[{i}] must be non-negative, got {dim}"
            )
        total *= dim
        if total > MAX_SAFE_SIZE:
            raise GenesisOverflowError(
                f"{name} total size ({total}) exceeds maximum ({MAX_SAFE_SIZE})"
            )

    return total


def validate_positive(value: int, name: str = "value") -> None:
    """Validate that a value is positive."""
    if value <= 0:
        raise GenesisValidationError(
            f"{name} must be positive, got {value}"
        )


def validate_non_negative(value: int, name: str = "value") -> None:
    """Validate that a value is non-negative."""
    if value < 0:
        raise GenesisValidationError(
            f"{name} must be non-negative, got {value}"
        )


def validate_pointer(ptr: ctypes.c_void_p, name: str = "pointer") -> None:
    """Validate that a pointer is not NULL."""
    if not ptr:
        raise GenesisMemoryError(f"{name} is NULL")


def validate_trajectory_dimensions(nframe: int, natom: int) -> None:
    """Validate trajectory dimensions won't overflow.

    Args:
        nframe: Number of frames
        natom: Number of atoms

    Raises:
        GenesisValidationError: If dimensions are invalid
        GenesisOverflowError: If total size exceeds safe limits
    """
    validate_non_negative(nframe, "nframe")
    validate_non_negative(natom, "natom")

    # Check nframe * natom * 3 * sizeof(double)
    total = nframe * natom * 3 * 8
    if total > MAX_SAFE_SIZE:
        raise GenesisOverflowError(
            f"Trajectory size ({nframe} frames × {natom} atoms × 3 × 8 = "
            f"{total} bytes) exceeds safe limit ({MAX_SAFE_SIZE})"
        )
