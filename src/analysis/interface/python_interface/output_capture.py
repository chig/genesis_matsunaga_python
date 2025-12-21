"""
Output capture utilities for GENESIS Python interface.

This module provides context managers for capturing stdout/stderr
from Fortran code while maintaining the ability to report errors.
"""

import os
import tempfile
from contextlib import contextmanager
from dataclasses import dataclass
from typing import Optional


@dataclass
class CapturedOutput:
    """Container for captured stdout and stderr content."""
    stdout: Optional[str] = None
    stderr: Optional[str] = None


@contextmanager
def capture_fortran_output(capture_stdout: bool = False,
                           capture_stderr: bool = True):
    """Context manager to capture Fortran output at file descriptor level.

    By default, suppresses stdout (to hide verbose Fortran output)
    and captures stderr (for error reporting).

    Args:
        capture_stdout: If True, capture stdout to a temp file.
                       If False, suppress stdout to /dev/null.
        capture_stderr: If True, capture stderr to a temp file.
                       If False, leave stderr unchanged.

    Yields:
        CapturedOutput: Object with captured stdout and stderr content
                       (populated after the context exits)

    Example:
        with capture_fortran_output() as output:
            # Call Fortran code
            ...
        if output.stderr:
            print(f"Fortran warnings: {output.stderr}")
    """
    saved_stdout = None
    saved_stderr = None
    stdout_file = None
    stderr_file = None
    result = CapturedOutput()

    try:
        # Handle stdout
        saved_stdout = os.dup(1)
        if capture_stdout:
            stdout_file = tempfile.NamedTemporaryFile(
                mode='w+', delete=False, suffix='.stdout'
            )
            os.dup2(stdout_file.fileno(), 1)
        else:
            # Suppress stdout to /dev/null
            devnull = os.open(os.devnull, os.O_WRONLY)
            os.dup2(devnull, 1)
            os.close(devnull)

        # Handle stderr
        if capture_stderr:
            saved_stderr = os.dup(2)
            stderr_file = tempfile.NamedTemporaryFile(
                mode='w+', delete=False, suffix='.stderr'
            )
            os.dup2(stderr_file.fileno(), 2)

        yield result

    finally:
        # Restore stdout
        if saved_stdout is not None:
            os.dup2(saved_stdout, 1)
            os.close(saved_stdout)

        # Restore stderr
        if saved_stderr is not None:
            os.dup2(saved_stderr, 2)
            os.close(saved_stderr)

        # Read captured stdout
        if stdout_file is not None:
            try:
                stdout_file.flush()
                stdout_file.seek(0)
                result.stdout = stdout_file.read()
            finally:
                stdout_file.close()
                try:
                    os.unlink(stdout_file.name)
                except OSError:
                    pass

        # Read captured stderr
        if stderr_file is not None:
            try:
                stderr_file.flush()
                stderr_file.seek(0)
                result.stderr = stderr_file.read()
            finally:
                stderr_file.close()
                try:
                    os.unlink(stderr_file.name)
                except OSError:
                    pass


@contextmanager
def suppress_stdout_capture_stderr():
    """Suppress stdout but capture stderr for error reporting.

    This is the recommended context manager for Fortran calls.
    It hides verbose Fortran output while preserving error messages.

    Yields:
        CapturedOutput: Object with captured stderr content

    Example:
        with suppress_stdout_capture_stderr() as output:
            LibGenesis().lib.some_analysis_c(...)
        if status.value != 0:
            raise GenesisFortranError(msg, stderr_output=output.stderr)
    """
    with capture_fortran_output(capture_stdout=False,
                                capture_stderr=True) as output:
        yield output


@contextmanager
def suppress_all_output():
    """Suppress both stdout and stderr.

    Use this when you want to completely silence Fortran output,
    including warnings. Note that this will also hide error messages
    that may be written to stderr.

    This is equivalent to the original suppress_stdout_simple().
    """
    saved_stdout = os.dup(1)
    saved_stderr = os.dup(2)

    try:
        devnull = os.open(os.devnull, os.O_WRONLY)
        os.dup2(devnull, 1)
        os.dup2(devnull, 2)
        os.close(devnull)

        yield

    finally:
        os.dup2(saved_stdout, 1)
        os.dup2(saved_stderr, 2)
        os.close(saved_stdout)
        os.close(saved_stderr)
