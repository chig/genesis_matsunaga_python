import os
from typing import IO
import numpy as np
import numpy.typing as npt


def read_msd_from_file(src_path: str | bytes | os.PathLike) \
        -> tuple(npt.NDArray[np.float64]):
    with open(src_path) as src:
        return read_msd_from_stream(src)


def read_msd_from_stream(src: IO[str]) -> tuple(npt.NDArray[np.float64]):
    lines = src.readlines()
    data_buf = []
    for line in lines:
        line_buf = []
        for v in line.split():
            line_buf.append(float(v))
        data_buf.append(line_buf)
    ndarray = np.array(data_buf)
    return ndarray

