#!/usr/bin/env python
"""
Unit tests for error handling in GENESIS Python interface.

This module tests the exception classes, validation utilities,
and output capture mechanisms.
"""

import ctypes
import unittest
import numpy as np

from .exceptions import (
    GenesisError,
    GenesisFortranError,
    GenesisValidationError,
    GenesisMemoryError,
    GenesisOverflowError,
)
from .validation import (
    validate_size,
    validate_pointer,
    validate_trajectory_dimensions,
    validate_positive,
    validate_non_negative,
    MAX_SAFE_SIZE,
)
from .output_capture import (
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
