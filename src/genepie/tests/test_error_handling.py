#!/usr/bin/env python
"""
Unit tests for error handling in GENESIS Python interface.

This module tests the exception classes, validation utilities,
and output capture mechanisms.
"""

import ctypes
import unittest
import numpy as np

from ..exceptions import (
    GenesisError,
    GenesisFortranError,
    GenesisFortranMemoryError,
    GenesisFortranFileError,
    GenesisFortranValidationError,
    GenesisFortranDataError,
    GenesisFortranNotSupportedError,
    GenesisFortranInternalError,
    GenesisValidationError,
    GenesisMemoryError,
    GenesisOverflowError,
    ErrorCode,
    raise_fortran_error,
)
from ..validation import (
    validate_size,
    validate_pointer,
    validate_trajectory_dimensions,
    validate_positive,
    validate_non_negative,
    MAX_SAFE_SIZE,
)
from ..output_capture import (
    capture_fortran_output,
    suppress_stdout_capture_stderr,
    suppress_all_output,
)


class TestExceptionHierarchy(unittest.TestCase):
    """Test exception class hierarchy."""

    def test_genesis_error_is_base(self):
        """All custom exceptions should inherit from GenesisError."""
        self.assertTrue(issubclass(GenesisFortranError, GenesisError))
        self.assertTrue(issubclass(GenesisValidationError, GenesisError))
        self.assertTrue(issubclass(GenesisMemoryError, GenesisError))
        self.assertTrue(issubclass(GenesisOverflowError, GenesisError))

    def test_genesis_error_inherits_from_exception(self):
        """GenesisError should inherit from Exception."""
        self.assertTrue(issubclass(GenesisError, Exception))

    def test_catch_all_genesis_errors(self):
        """Should be able to catch all errors with GenesisError."""
        errors = [
            GenesisFortranError("test"),
            GenesisValidationError("test"),
            GenesisMemoryError("test"),
            GenesisOverflowError("test"),
        ]
        for error in errors:
            with self.assertRaises(GenesisError):
                raise error


class TestErrorCode(unittest.TestCase):
    """Test ErrorCode class."""

    def test_memory_error_category(self):
        """Memory errors (100-199) should return MEMORY_ERROR."""
        self.assertEqual(ErrorCode.get_category(101), "MEMORY_ERROR")
        self.assertEqual(ErrorCode.get_category(102), "MEMORY_ERROR")
        self.assertEqual(ErrorCode.get_category(199), "MEMORY_ERROR")

    def test_file_error_category(self):
        """File errors (200-299) should return FILE_ERROR."""
        self.assertEqual(ErrorCode.get_category(201), "FILE_ERROR")
        self.assertEqual(ErrorCode.get_category(299), "FILE_ERROR")

    def test_validation_error_category(self):
        """Validation errors (300-399) should return VALIDATION_ERROR."""
        self.assertEqual(ErrorCode.get_category(301), "VALIDATION_ERROR")
        self.assertEqual(ErrorCode.get_category(399), "VALIDATION_ERROR")

    def test_data_error_category(self):
        """Data errors (400-499) should return DATA_ERROR."""
        self.assertEqual(ErrorCode.get_category(401), "DATA_ERROR")
        self.assertEqual(ErrorCode.get_category(499), "DATA_ERROR")

    def test_not_supported_error_category(self):
        """Not supported errors (500-599) should return NOT_SUPPORTED_ERROR."""
        self.assertEqual(ErrorCode.get_category(501), "NOT_SUPPORTED_ERROR")
        self.assertEqual(ErrorCode.get_category(599), "NOT_SUPPORTED_ERROR")

    def test_internal_error_category(self):
        """Internal errors (600-699) should return INTERNAL_ERROR."""
        self.assertEqual(ErrorCode.get_category(601), "INTERNAL_ERROR")
        self.assertEqual(ErrorCode.get_category(699), "INTERNAL_ERROR")

    def test_unknown_error_category(self):
        """Unknown codes should return UNKNOWN_ERROR."""
        self.assertEqual(ErrorCode.get_category(0), "UNKNOWN_ERROR")
        self.assertEqual(ErrorCode.get_category(99), "UNKNOWN_ERROR")
        self.assertEqual(ErrorCode.get_category(700), "UNKNOWN_ERROR")

    def test_error_code_constants(self):
        """Error code constants should be defined."""
        self.assertEqual(ErrorCode.ERROR_ALLOC, 101)
        self.assertEqual(ErrorCode.ERROR_FILE_NOT_FOUND, 201)
        self.assertEqual(ErrorCode.ERROR_INVALID_PARAM, 301)
        self.assertEqual(ErrorCode.ERROR_DATA_MISMATCH, 401)
        self.assertEqual(ErrorCode.ERROR_NOT_SUPPORTED, 501)
        self.assertEqual(ErrorCode.ERROR_INTERNAL, 601)


class TestFortranErrorSubclasses(unittest.TestCase):
    """Test Fortran error subclasses."""

    def test_memory_error_inheritance(self):
        """GenesisFortranMemoryError should inherit from GenesisFortranError."""
        self.assertTrue(issubclass(GenesisFortranMemoryError, GenesisFortranError))
        error = GenesisFortranMemoryError("alloc failed", code=101)
        self.assertIsInstance(error, GenesisFortranError)
        self.assertIsInstance(error, GenesisError)

    def test_file_error_inheritance(self):
        """GenesisFortranFileError should inherit from GenesisFortranError."""
        self.assertTrue(issubclass(GenesisFortranFileError, GenesisFortranError))

    def test_validation_error_inheritance(self):
        """GenesisFortranValidationError should inherit from GenesisFortranError."""
        self.assertTrue(issubclass(GenesisFortranValidationError, GenesisFortranError))

    def test_data_error_inheritance(self):
        """GenesisFortranDataError should inherit from GenesisFortranError."""
        self.assertTrue(issubclass(GenesisFortranDataError, GenesisFortranError))

    def test_not_supported_error_inheritance(self):
        """GenesisFortranNotSupportedError should inherit from GenesisFortranError."""
        self.assertTrue(issubclass(GenesisFortranNotSupportedError, GenesisFortranError))

    def test_internal_error_inheritance(self):
        """GenesisFortranInternalError should inherit from GenesisFortranError."""
        self.assertTrue(issubclass(GenesisFortranInternalError, GenesisFortranError))

    def test_catch_all_with_base_class(self):
        """All subclasses should be catchable with GenesisFortranError."""
        errors = [
            GenesisFortranMemoryError("test", code=101),
            GenesisFortranFileError("test", code=201),
            GenesisFortranValidationError("test", code=301),
            GenesisFortranDataError("test", code=401),
            GenesisFortranNotSupportedError("test", code=501),
            GenesisFortranInternalError("test", code=601),
        ]
        for error in errors:
            with self.assertRaises(GenesisFortranError):
                raise error


class TestRaiseFortranError(unittest.TestCase):
    """Test raise_fortran_error factory function."""

    def test_raises_memory_error(self):
        """Code 101 should raise GenesisFortranMemoryError."""
        with self.assertRaises(GenesisFortranMemoryError) as ctx:
            raise_fortran_error("test", code=101)
        self.assertEqual(ctx.exception.code, 101)

    def test_raises_file_error(self):
        """Code 201 should raise GenesisFortranFileError."""
        with self.assertRaises(GenesisFortranFileError) as ctx:
            raise_fortran_error("test", code=201)
        self.assertEqual(ctx.exception.code, 201)

    def test_raises_validation_error(self):
        """Code 301 should raise GenesisFortranValidationError."""
        with self.assertRaises(GenesisFortranValidationError) as ctx:
            raise_fortran_error("test", code=301)
        self.assertEqual(ctx.exception.code, 301)

    def test_raises_data_error(self):
        """Code 401 should raise GenesisFortranDataError."""
        with self.assertRaises(GenesisFortranDataError) as ctx:
            raise_fortran_error("test", code=401)
        self.assertEqual(ctx.exception.code, 401)

    def test_raises_not_supported_error(self):
        """Code 501 should raise GenesisFortranNotSupportedError."""
        with self.assertRaises(GenesisFortranNotSupportedError) as ctx:
            raise_fortran_error("test", code=501)
        self.assertEqual(ctx.exception.code, 501)

    def test_raises_internal_error(self):
        """Code 601 should raise GenesisFortranInternalError."""
        with self.assertRaises(GenesisFortranInternalError) as ctx:
            raise_fortran_error("test", code=601)
        self.assertEqual(ctx.exception.code, 601)

    def test_raises_base_for_unknown_code(self):
        """Unknown code should raise GenesisFortranError."""
        with self.assertRaises(GenesisFortranError) as ctx:
            raise_fortran_error("test", code=999)
        self.assertEqual(ctx.exception.code, 999)
        # Should not be a subclass
        self.assertEqual(type(ctx.exception), GenesisFortranError)

    def test_preserves_message_and_stderr(self):
        """Should preserve message and stderr_output."""
        with self.assertRaises(GenesisFortranError) as ctx:
            raise_fortran_error("error message", code=101, stderr_output="stderr text")
        self.assertEqual(str(ctx.exception), "error message")
        self.assertEqual(ctx.exception.stderr_output, "stderr text")


class TestGenesisFortranErrorCategoryMethods(unittest.TestCase):
    """Test is_*_error() methods on GenesisFortranError."""

    def test_is_memory_error(self):
        """is_memory_error() should return True for code 100-199."""
        error = GenesisFortranError("test", code=101)
        self.assertTrue(error.is_memory_error())
        self.assertFalse(error.is_file_error())

    def test_is_file_error(self):
        """is_file_error() should return True for code 200-299."""
        error = GenesisFortranError("test", code=201)
        self.assertTrue(error.is_file_error())
        self.assertFalse(error.is_memory_error())

    def test_is_not_supported_error(self):
        """is_not_supported_error() should return True for code 500-599."""
        error = GenesisFortranError("test", code=502)
        self.assertTrue(error.is_not_supported_error())
        self.assertFalse(error.is_internal_error())

    def test_category_attribute(self):
        """category attribute should be set correctly."""
        error = GenesisFortranError("test", code=501)
        self.assertEqual(error.category, "NOT_SUPPORTED_ERROR")


class TestGenesisFortranError(unittest.TestCase):
    """Test GenesisFortranError exception."""

    def test_basic_creation(self):
        """Test creating error with just a message."""
        error = GenesisFortranError("Test error")
        self.assertEqual(str(error), "Test error")
        self.assertEqual(error.code, 0)
        self.assertEqual(error.stderr_output, "")

    def test_with_code(self):
        """Test creating error with code."""
        error = GenesisFortranError("Test error", code=101)
        self.assertEqual(error.code, 101)

    def test_with_stderr(self):
        """Test creating error with stderr output."""
        error = GenesisFortranError(
            "Test error",
            code=101,
            stderr_output="Warning: something happened"
        )
        self.assertEqual(error.stderr_output, "Warning: something happened")

    def test_empty_stderr(self):
        """Empty stderr should be stored as empty string."""
        error = GenesisFortranError("Test error", stderr_output="")
        self.assertEqual(error.stderr_output, "")

    def test_none_stderr(self):
        """None stderr should be converted to empty string."""
        error = GenesisFortranError("Test error", stderr_output=None)
        self.assertEqual(error.stderr_output, "")


class TestValidateSize(unittest.TestCase):
    """Test validate_size function."""

    def test_valid_positive_int(self):
        """Valid positive integer should return the value."""
        result = validate_size(100, "test")
        self.assertEqual(result, 100)

    def test_valid_zero(self):
        """Zero is a valid size."""
        result = validate_size(0, "test")
        self.assertEqual(result, 0)

    def test_negative_raises_error(self):
        """Negative size should raise GenesisValidationError."""
        with self.assertRaises(GenesisValidationError) as ctx:
            validate_size(-1, "array_size")
        self.assertIn("array_size", str(ctx.exception))
        self.assertIn("non-negative", str(ctx.exception))

    def test_valid_sequence(self):
        """Valid sequence should return product."""
        result = validate_size([10, 20, 30], "test")
        self.assertEqual(result, 6000)

    def test_sequence_with_negative(self):
        """Sequence with negative dimension should raise error."""
        with self.assertRaises(GenesisValidationError):
            validate_size([10, -5, 30], "test")

    def test_overflow_detection(self):
        """Should detect overflow for very large sizes."""
        with self.assertRaises(GenesisOverflowError):
            validate_size(MAX_SAFE_SIZE + 1, "test")

    def test_sequence_overflow_detection(self):
        """Should detect overflow in sequence products."""
        # This product would overflow
        with self.assertRaises(GenesisOverflowError):
            validate_size([100000, 100000, 100000], "test")


class TestValidatePointer(unittest.TestCase):
    """Test validate_pointer function."""

    def test_null_pointer_raises_error(self):
        """NULL pointer should raise GenesisMemoryError."""
        with self.assertRaises(GenesisMemoryError) as ctx:
            validate_pointer(None, "my_pointer")
        self.assertIn("my_pointer", str(ctx.exception))

    def test_valid_pointer_passes(self):
        """Valid pointer should not raise."""
        ptr = ctypes.c_void_p(1234)
        validate_pointer(ptr, "test")  # Should not raise


class TestValidateTrajectoryDimensions(unittest.TestCase):
    """Test validate_trajectory_dimensions function."""

    def test_valid_dimensions(self):
        """Valid dimensions should pass."""
        validate_trajectory_dimensions(1000, 10000)  # 30M elements, OK

    def test_zero_dimensions(self):
        """Zero dimensions should pass."""
        validate_trajectory_dimensions(0, 100)
        validate_trajectory_dimensions(100, 0)

    def test_negative_nframe(self):
        """Negative nframe should raise error."""
        with self.assertRaises(GenesisValidationError):
            validate_trajectory_dimensions(-1, 100)

    def test_negative_natom(self):
        """Negative natom should raise error."""
        with self.assertRaises(GenesisValidationError):
            validate_trajectory_dimensions(100, -1)

    def test_overflow_detection(self):
        """Should detect overflow for very large trajectories."""
        with self.assertRaises(GenesisOverflowError):
            validate_trajectory_dimensions(100000, 100000)


class TestValidatePositive(unittest.TestCase):
    """Test validate_positive function."""

    def test_positive_passes(self):
        """Positive value should pass."""
        validate_positive(1, "test")
        validate_positive(100, "test")

    def test_zero_fails(self):
        """Zero should raise error."""
        with self.assertRaises(GenesisValidationError):
            validate_positive(0, "test")

    def test_negative_fails(self):
        """Negative should raise error."""
        with self.assertRaises(GenesisValidationError):
            validate_positive(-1, "test")


class TestValidateNonNegative(unittest.TestCase):
    """Test validate_non_negative function."""

    def test_positive_passes(self):
        """Positive value should pass."""
        validate_non_negative(1, "test")
        validate_non_negative(100, "test")

    def test_zero_passes(self):
        """Zero should pass."""
        validate_non_negative(0, "test")

    def test_negative_fails(self):
        """Negative should raise error."""
        with self.assertRaises(GenesisValidationError):
            validate_non_negative(-1, "test")


class TestOutputCapture(unittest.TestCase):
    """Test output capture context managers."""

    def test_suppress_stdout_capture_stderr_basic(self):
        """Basic test of suppress_stdout_capture_stderr."""
        with suppress_stdout_capture_stderr() as captured:
            pass  # No output
        # stdout is suppressed, not captured, so it should be None
        self.assertIsNone(captured.stdout)
        # stderr is captured (may be empty string)
        self.assertIsInstance(captured.stderr, (str, type(None)))

    def test_suppress_all_output(self):
        """Test suppress_all_output suppresses everything."""
        with suppress_all_output():
            print("This should be suppressed")
            # No error should occur


if __name__ == "__main__":
    unittest.main()
