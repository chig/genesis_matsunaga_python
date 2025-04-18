import ctypes


class STrajectoriesC(ctypes.Structure):
    """
    Definition for accessing the Genesis's s_trajectories object
    allocated in the memory area of the C language interface from Python.
    """
    _fields_ = [("coords", ctypes.c_void_p),
                ("pbc_boxes", ctypes.c_void_p),
                ("nframe", ctypes.c_int),
                ("natom", ctypes.c_int)]
