import ctypes
import numpy as np
import numpy.typing as npt


def write_bool_ndarray(src: npt.NDArray[np.bool_], dst: ctypes.c_void_p):
    dst_p = ctypes.cast(dst, ctypes.POINTER(ctypes.c_bool))
    for i, v in enumerate(np.ravel(src)):
        dst_p[i] = ctypes.c_bool(v)


def write_int_ndarray(src: npt.NDArray[np.int64], dst: ctypes.c_void_p):
    dst_p = ctypes.cast(dst, ctypes.POINTER(ctypes.c_int))
    for i, v in enumerate(np.ravel(src)):
        dst_p[i] = ctypes.c_int(v)


def write_double_ndarray(src: npt.NDArray[np.float64], dst: ctypes.c_void_p):
    dst_p = ctypes.cast(dst, ctypes.POINTER(ctypes.c_double))
    for i, v in enumerate(np.ravel(src)):
        dst_p[i] = ctypes.c_double(v)


def write_fixed_length_string_ndarray(src: npt.NDArray[np.str_], dst: ctypes.c_void_p):
    dst_p = ctypes.cast(dst, ctypes.POINTER(ctypes.c_char))
    for i, v in enumerate(np.ravel(src)):
        dst_p[i] = v.encode('ascii')
