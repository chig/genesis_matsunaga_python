import ctypes
from typing import Sequence, Union

import numpy as np
import numpy.typing as npt

from .exceptions import GenesisValidationError, GenesisMemoryError
from .validation import validate_size, validate_pointer


def conv_ndarray(src: ctypes.c_void_p, size: Union[int, Sequence[int]],
                 type_clang: type, type_numpy: type) -> npt.NDArray:
    if not src:
        return np.empty(0, dtype=type_numpy)

    # Validate size before iteration
    all_size = validate_size(size, "array size")
    if all_size == 0:
        return np.empty(0, dtype=type_numpy)

    # Normalize size to tuple
    if isinstance(size, (int, np.integer)):
        size = (int(size),)
    else:
        size = tuple(size)

    ptr = ctypes.cast(src, ctypes.POINTER(type_clang))

    # Use count parameter to bound iteration
    return (np.fromiter((ptr[i] for i in range(all_size)),
                        dtype=type_numpy, count=all_size)
            .reshape(size))


def conv_bool_ndarray(src: ctypes.c_void_p, size: Union[int, Sequence[int]]) \
        -> npt.NDArray[np.bool_]:
    """Convert C bool array to numpy bool array."""
    return conv_ndarray(src, size, ctypes.c_bool, np.bool_)


def conv_int_ndarray(src: ctypes.c_void_p, size: Union[int, Sequence[int]]) \
        -> npt.NDArray[np.int64]:
    """Convert C int array to numpy int64 array."""
    return conv_ndarray(src, size, ctypes.c_int, np.int64)


def conv_double_ndarray(src: ctypes.c_void_p, size: Union[int, Sequence[int]]) \
        -> npt.NDArray[np.float64]:
    """Convert C double array to numpy float64 array."""
    return conv_ndarray(src, size, ctypes.c_double, np.float64)


def conv_fixed_length_string_ndarray(
        str_array: ctypes.c_void_p, size: Union[int, Sequence[int]]) \
                -> npt.NDArray[np.str_]:
    """Convert C fixed-length string array to numpy string array."""
    if not str_array:
        return np.empty(0, dtype=np.str_)

    # Validate size before iteration
    all_size = validate_size(size, "string array size")
    if all_size == 0:
        return np.empty(0, dtype=np.str_)

    # Normalize size to tuple
    if isinstance(size, (int, np.integer)):
        size = (int(size),)
    else:
        size = tuple(size)

    ptr = ctypes.cast(str_array, ctypes.POINTER(ctypes.c_char))
    dst = np.empty(all_size, dtype=str)

    for i in range(all_size):
        try:
            dst[i] = ptr[i].decode('ascii')
        except (UnicodeDecodeError, ValueError) as e:
            raise GenesisValidationError(
                f"Failed to decode character at index {i}: {e}"
            )

    return dst.reshape(size)


def conv_pystring_ndarray(
        str_array: ctypes.c_void_p, size: tuple[int, int]) \
                -> npt.NDArray[np.object_]:
    """Convert C string array to numpy object array of Python strings."""
    if not str_array:
        return np.empty(0, dtype=np.object_)

    # Validate dimensions
    if len(size) != 2:
        raise GenesisValidationError(
            f"Expected 2D size tuple, got {len(size)}D"
        )

    num_strings, str_len = size
    validate_size(num_strings, "number of strings")
    validate_size(str_len, "string length")

    if num_strings == 0:
        return np.empty(0, dtype=np.object_)

    # Validate total size for pointer arithmetic
    validate_size(num_strings * str_len, "total buffer size")

    ptr = ctypes.cast(str_array, ctypes.POINTER(ctypes.c_char))
    dst = np.empty(num_strings, dtype=np.object_)

    for i in range(num_strings):
        offset = i * str_len
        cur_ptr = ctypes.cast(
            ctypes.addressof(ptr.contents) + offset,
            ctypes.POINTER(ctypes.c_char)
        )
        dst[i] = conv_string(cur_ptr, str_len)

    return dst


def conv_string(c_char_ptr: ctypes.c_void_p, size: int = -1,
                encoding: str = 'utf-8') -> str:
    """Convert C string to Python string."""
    if not c_char_ptr:
        raise GenesisMemoryError("c_char_ptr is NULL")

    if size > 0:
        validate_size(size, "string size")

    byte_str = ctypes.string_at(c_char_ptr, size)
    return byte_str.decode(encoding, errors='replace').rstrip('\x00')
