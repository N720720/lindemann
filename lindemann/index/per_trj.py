from typing import Any

import bottleneck as bn
import numba as nb
import numpy as np
import numpy.typing as npt

# No typing for jit functions https://github.com/numba/numba/issues/7424


@nb.njit(fastmath=True, error_model="numpy")  # type: ignore
def lindemann_per_atom(frames: npt.NDArray[np.float32]) -> Any:

    """Calculate the lindeman index
    Args:
        frames: numpy array of shape(frames,atoms)
    Returns:
        float32: returns the lindeman index
    """

    dt = frames.dtype
    natoms = len(frames[0])
    nframes = len(frames)
    array_mean = np.zeros((natoms, natoms), dtype=dt)
    array_var = np.zeros((natoms, natoms), dtype=dt)
    iframe = dt.type(1)
    for coords in frames:

        # here we do someting similar to scipy's spatial.distance.pdist scipy.spatial.distance.pdist
        n, p = coords.shape
        array_distance = np.zeros((n, n), dtype=dt)
        for i in range(n):
            for j in range(i + 1, n):
                d = dt.type(0.0)
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
                array_mean[i, j] = mean + delta / iframe
                array_var[i, j] = var + delta * (xn - array_mean[i, j])
        iframe += 1.0  # type: ignore[assignment]
        if iframe > nframes:
            break

    for i in range(natoms):
        for j in range(i + 1, natoms):
            array_mean[j, i] = array_mean[i, j]
            array_var[j, i] = array_var[i, j]

    lindemann_indices = np.divide(np.sqrt(np.divide(array_var, nframes)), array_mean)
    return lindemann_indices


def calculate(frames: npt.NDArray[np.float64]) -> float:

    return np.mean(bn.nanmean(lindemann_per_atom(frames), axis=1))  # type: ignore[no-any-return, no-untyped-call]
