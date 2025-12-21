"""
Custom exception classes for GENESIS Python interface.
"""

from typing import Optional


class GenesisError(Exception):
    """Base exception for all GENESIS Python interface errors."""
    pass


class GenesisValidationError(GenesisError):
    """Raised when input validation fails on the Python side."""
    pass


class GenesisFortranError(GenesisError):
    """Raised when Fortran code returns an error status."""

    def __init__(self, message: str, code: int = 0,
                 stderr_output: Optional[str] = None):
        self.code = code
        self.stderr_output = stderr_output or ""
        super().__init__(message)


class GenesisMemoryError(GenesisError):
    """Raised when memory allocation or pointer operations fail."""
    pass


class GenesisOverflowError(GenesisError):
    """Raised when size calculations would overflow."""
    pass
