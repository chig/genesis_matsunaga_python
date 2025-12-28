"""
Custom exception classes for GENESIS Python interface.
"""

from typing import Optional


class ErrorCode:
    """Error code constants matching Fortran error_mod.fpp."""

    # Memory errors (100-199)
    ERROR_ALLOC = 101
    ERROR_DEALLOC = 102
    ERROR_NOT_ALLOCATED = 103

    # File errors (200-299)
    ERROR_FILE_NOT_FOUND = 201
    ERROR_FILE_FORMAT = 202
    ERROR_FILE_READ = 203

    # Validation errors (300-399)
    ERROR_INVALID_PARAM = 301
    ERROR_MISSING_PARAM = 302
    ERROR_DIMENSION = 303
    ERROR_ATOM_COUNT = 304
    ERROR_GRID_SIZE = 305
    ERROR_BOND_INFO = 306
    ERROR_SELECTION = 307

    # Data errors (400-499)
    ERROR_DATA_MISMATCH = 401
    ERROR_NO_DATA = 402
    ERROR_MASS_UNDEFINED = 403
    ERROR_PBC_BOX = 404

    # Not supported errors (500-599)
    ERROR_NOT_SUPPORTED = 501
    ERROR_DIM_NOT_SUPP = 502
    ERROR_FUNC_NOT_SUPP = 503
    ERROR_BLOCK_NOT_SUPP = 504

    # Internal errors (600-699)
    ERROR_INTERNAL = 601
    ERROR_SYNTAX = 602
    ERROR_GENERIC = 699

    @classmethod
    def get_category(cls, code: int) -> str:
        """Get the error category name from a code."""
        if 100 <= code < 200:
            return "MEMORY_ERROR"
        elif 200 <= code < 300:
            return "FILE_ERROR"
        elif 300 <= code < 400:
            return "VALIDATION_ERROR"
        elif 400 <= code < 500:
            return "DATA_ERROR"
        elif 500 <= code < 600:
            return "NOT_SUPPORTED_ERROR"
        elif 600 <= code < 700:
            return "INTERNAL_ERROR"
        else:
            return "UNKNOWN_ERROR"


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
        self.category = ErrorCode.get_category(code)
        super().__init__(message)

    def is_memory_error(self) -> bool:
        """Check if this is a memory-related error."""
        return 100 <= self.code < 200

    def is_file_error(self) -> bool:
        """Check if this is a file-related error."""
        return 200 <= self.code < 300

    def is_validation_error(self) -> bool:
        """Check if this is a validation error."""
        return 300 <= self.code < 400

    def is_data_error(self) -> bool:
        """Check if this is a data-related error."""
        return 400 <= self.code < 500

    def is_not_supported_error(self) -> bool:
        """Check if this is a not-supported error."""
        return 500 <= self.code < 600

    def is_internal_error(self) -> bool:
        """Check if this is an internal error."""
        return 600 <= self.code < 700


# Specific Fortran error subclasses for type-based exception handling
class GenesisFortranMemoryError(GenesisFortranError):
    """Fortran memory allocation/deallocation error."""
    pass


class GenesisFortranFileError(GenesisFortranError):
    """Fortran file I/O error."""
    pass


class GenesisFortranValidationError(GenesisFortranError):
    """Fortran input validation error."""
    pass


class GenesisFortranDataError(GenesisFortranError):
    """Fortran data mismatch error."""
    pass


class GenesisFortranNotSupportedError(GenesisFortranError):
    """Fortran feature not supported error."""
    pass


class GenesisFortranInternalError(GenesisFortranError):
    """Fortran internal error."""
    pass


class GenesisMemoryError(GenesisError):
    """Raised when memory allocation or pointer operations fail."""
    pass


class GenesisOverflowError(GenesisError):
    """Raised when size calculations would overflow."""
    pass


def raise_fortran_error(
    message: str,
    code: int = 0,
    stderr_output: Optional[str] = None
) -> None:
    """
    Factory function to raise the appropriate exception subclass.

    This allows Python code to catch specific error types:

        try:
            result = genesis_exe.rmsd_analysis(...)
        except GenesisFortranNotSupportedError as e:
            print(f"Feature not supported: {e}")
        except GenesisFortranValidationError as e:
            print(f"Invalid input: {e}")
        except GenesisFortranError as e:
            print(f"Fortran error: {e}")

    Args:
        message: Error message
        code: Fortran error code
        stderr_output: Captured stderr output from Fortran
    """
    if 100 <= code < 200:
        raise GenesisFortranMemoryError(message, code, stderr_output)
    elif 200 <= code < 300:
        raise GenesisFortranFileError(message, code, stderr_output)
    elif 300 <= code < 400:
        raise GenesisFortranValidationError(message, code, stderr_output)
    elif 400 <= code < 500:
        raise GenesisFortranDataError(message, code, stderr_output)
    elif 500 <= code < 600:
        raise GenesisFortranNotSupportedError(message, code, stderr_output)
    elif 600 <= code < 700:
        raise GenesisFortranInternalError(message, code, stderr_output)
    else:
        raise GenesisFortranError(message, code, stderr_output)
