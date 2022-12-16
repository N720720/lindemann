from typing import List

import numba as nb
import numpy as np
import numpy.typing as npt
from numba import float32


@nb.njit(fastmath=True, error_model="numpy")  # type: ignore # , cache=True) #(parallel=True)
def calculate(frames: npt.NDArray[np.float32]) -> npt.NDArray[np.float32]:

    """Calculate the contribution of each atom to the lindemann index over the frames

    Args:
        frames: numpy array of shape(frames,atoms)
    Returns:
        npt.NDArray[np.float32]: Returns 1D array with the progression of the lindeman index per frame of shape(frames, atoms)
    """

    first = True
    # natoms = natoms
    dt = frames.dtype
    natoms = len(frames[0])
    nframes = len(frames)
    len_frames = len(frames)
    array_mean = np.zeros((natoms, natoms), dtype=dt)
    array_var = np.zeros((natoms, natoms), dtype=dt)
    iframe = dt.type(1)
    lindex_array = np.zeros((len_frames, natoms), dtype=dt)
    for q, coords in enumerate(frames):
        n, p = coords.shape
        array_distance = np.zeros((n, n), dtype=dt)
        for i in range(n):
            for j in range(i + 1, n):
                d = 0.0
                for k in range(p):
                    d += (coords[i, k] - coords[j, k]) ** dt.type(2)
                array_distance[i, j] = np.sqrt(d)
        array_distance += array_distance.T

        #################################################################################
        # update mean and var arrays based on Welford algorithm suggested by Donald Knuth
        #################################################################################
        for i in range(natoms):
            for j in range(i + 1, natoms):
                xn = array_distance[i, j]
                mean = array_mean[i, j]
                var = array_var[i, j]
                delta = xn - mean
                # update mean
                array_mean[i, j] = mean + delta / iframe
                # update variance
                array_var[i, j] = var + delta * (xn - array_mean[i, j])
        iframe += 1  # type: ignore[assignment]
        if iframe > nframes + 1:
            break

        for i in range(natoms):
            for j in range(i + 1, natoms):
                array_mean[j, i] = array_mean[i, j]
                array_var[j, i] = array_var[i, j]

        if first:
            lindemann_indices = np.zeros((natoms), dtype=dt)
            first = False
        else:
            np.fill_diagonal(array_mean, 1)
            lindemann_indices = np.zeros((natoms), dtype=dt)
            lindemann_indices = np.divide(np.sqrt(np.divide(array_var, iframe - 1)), array_mean)
            lindemann_indices = np.asarray([np.mean(lin[lin != 0]) for lin in lindemann_indices])

        lindex_array[q] = lindemann_indices
    return lindex_array