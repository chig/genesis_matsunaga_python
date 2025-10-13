# -*- coding: utf-8 -*-
"""
genesis_log.py
A small parser for GENESIS MD log files (e.g., log1).
- Parses "INFO" header lines to build a column index (STEP, TEMPERATURE, ...).
- Reads following numeric lines using that index.

Compatible with Python 3.8+ (tested up to 3.13).
"""

from __future__ import annotations
from pathlib import Path
import re
from typing import Dict, Iterable, List, Tuple, Union

_INFO_RE = re.compile(r"^INFO\.")  # lines starting with "INFO."
_WS_NL_RE = re.compile(r"\n$")

def _iter_lines(fp: Iterable[str]):
    """Yield stripped lines as-is (no trailing newline)."""
    for line in fp:
        yield _WS_NL_RE.sub("", line)

def _is_header_tokens(tokens: List[str]) -> bool:
    """
    Heuristic used in the original script: treat as header if any token ends with an alphabetic char.
    """
    return any(tok and tok[-1].isalpha() for tok in tokens)

def _build_index_from_header(tokens: List[str]) -> Dict[str, int]:
    """
    Map header tokens to their column indices (STEP, TEMPERATURE, PRESSURE, ...).
    """
    index: Dict[str, int] = {}
    for i, tok in enumerate(tokens):
        if tok and tok[-1].isalpha():
            index[tok] = i
    return index

def _safe_float(x: str) -> float:
    # be tolerant of things like "1.23E+03", "nan", etc.
    return float(x)

def read(
    filename: Union[str, Path],
    keyword: str = "TEMPERATURE",
    step_key: str = "STEP",
) -> Tuple[List[float], List[float]]:
    """
    Read (steps, values) for a given keyword from a GENESIS log.

    Parameters
    ----------
    filename : str or Path
        Log file path (e.g., "log1").
    keyword : str
        Column to extract (e.g., "TEMPERATURE", "ENERGY", "PRESSURE", ...).
    step_key : str
        Column used as x-axis (default: "STEP").

    Returns
    -------
    (x, y) : (List[float], List[float])
        Steps and the corresponding values for `keyword`.
    """
    path = Path(filename)
    if not path.exists():
        raise FileNotFoundError(f"GENESIS log not found: {path}")

    x: List[float] = []
    y: List[float] = []
    index: Dict[str, int] = {}

    with path.open("r", encoding="utf-8", errors="replace") as f:
        for raw in _iter_lines(f):
            if not raw:
                continue

            # Only parse lines starting with INFO.
            if not _INFO_RE.match(raw):
                continue

            # Drop "INFO."
            payload = _INFO_RE.sub("", raw).strip()
            tokens = payload.split()

            if not tokens:
                continue

            if _is_header_tokens(tokens):
                # Header line -> (re)build index
                index = _build_index_from_header(tokens)
                continue

            # Data line -> use current index (if available)
            if not index:
                # No header yet; skip to be safe
                continue

            if step_key not in index or keyword not in index:
                # Column not present in current header
                continue

            try:
                step = _safe_float(tokens[index[step_key]])
                val  = _safe_float(tokens[index[keyword]])
            except (IndexError, ValueError):
                # Malformed row; skip
                continue

            x.append(step)
            y.append(val)

    return x, y

# ---- Optional utilities -----------------------------------------------------

def read_many(
    filename: Union[str, Path],
    keys: Iterable[str],
    step_key: str = "STEP",
) -> Tuple[List[float], Dict[str, List[float]]]:
    """
    Read multiple columns at once. Returns (steps, {key: values}).
    """
    steps: List[float] = []
    series: Dict[str, List[float]] = {k: [] for k in keys}

    # Reuse the single-key reader to avoid code duplication
    # First pass to get steps from step_key (or TEMPERATURE as a proxy if step only appears there)
    # Here we just parse once and fill all keys together.
    path = Path(filename)
    if not path.exists():
        raise FileNotFoundError(f"GENESIS log not found: {path}")

    index: Dict[str, int] = {}
    with path.open("r", encoding="utf-8", errors="replace") as f:
        for raw in _iter_lines(f):
            if not raw:
                continue
            if not _INFO_RE.match(raw):
                continue

            payload = _INFO_RE.sub("", raw).strip()
            tokens = payload.split()
            if not tokens:
                continue

            if _is_header_tokens(tokens):
                index = _build_index_from_header(tokens)
                continue

            if not index or (step_key not in index):
                continue

            try:
                step = _safe_float(tokens[index[step_key]])
            except (IndexError, ValueError):
                continue

            # Collect for all requested keys that exist in the header
            row_vals = {}
            any_key_ok = False
            for k in keys:
                if k in index:
                    try:
                        row_vals[k] = _safe_float(tokens[index[k]])
                        any_key_ok = True
                    except (IndexError, ValueError):
                        row_vals[k] = None
                else:
                    row_vals[k] = None

            if not any_key_ok:
                continue

            steps.append(step)
            for k in keys:
                series[k].append(row_vals[k])

    return steps, series

# ---- Backward-compatible wrapper -------------------------------------------

def read_genesis(filename: str, key: str):
    """
    Back-compat wrapper to keep the tutorial call signature:
        (x, y) = read_genesis("log1", "TEMPERATURE")
    """
    return read(filename, keyword=key)

