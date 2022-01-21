from typing import Any

import numba as nb
import numpy as np
import numpy.typing as npt

# No typing for jit functions https://github.com/numba/numba/issues/7424


@nb.njit
def np_apply_along_axis(func1d, axis, arr):
    assert arr.ndim == 2
    assert axis in [0, 1]
    if axis == 0:
        result = np.empty(arr.shape[1])
        for i in range(len(result)):
            result[i] = func1d(arr[:, i])
    else:
        result = np.empty(arr.shape[0])
        for i in range(len(result)):
            result[i] = func1d(arr[i, :])
    return result


@nb.njit
def np_mean(array, axis):
    return np_apply_along_axis(np.mean, axis, array)


@nb.njit
def np_std(array, axis):
    return np_apply_along_axis(np.std, axis, array)


@nb.njit(fastmath=True, error_model="numpy")  # type: ignore
def calculate(frames: npt.NDArray[np.float64]) -> Any:

    """Calculates the lindemann index for """
    dt = frames.dtype
    natoms = len(frames[0])
    nframes = len(frames)
    array_mean = np.zeros((natoms, natoms), dtype=dt)
    array_var = np.zeros((natoms, natoms), dtype=dt)
    # array_distance = np.zeros((natoms, natoms),dtype=dt)
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
        iframe += 1.0
        if iframe > nframes:
            break

    for i in range(natoms):
        for j in range(i + 1, natoms):
            array_mean[j, i] = array_mean[i, j]
            array_var[j, i] = array_var[i, j]
    np.fill_diagonal(array_mean, 1)
    # array_mean = np.array([array_mean[i][array_mean[i] != array_mean[i][i]] for i in range(len(array_mean))])
    # array_mean = np.array([array_var[i][array_var[i] != array_var[i][i]] for i in range(len(array_var))])
    lindemann_indices = np.divide(np.sqrt(np.divide(array_var, nframes)), array_mean)
    # np.fill_diagonal(lindemann_indices,0)
    # zap = np.zeros((459), dtype=dt)
    # zap =
    # zappo = np.mean(np.asarray(zap))

    return np.mean(
        np.asarray([np.mean(lin[lin != 0]) for lin in lindemann_indices])
    )  # , a#a#np.mean(np_mean(lindemann_indices,axis=1))


# def calculate(frames: npt.NDArray[np.float64]) -> float:
#     """
#     Small helper function, since numba has not implemented the np.nanmean with axis parameter
#     I cant implemnet this in the jit function for now.
#     """
#     zap = lindemann_per_atom(frames)
#     return zap#np.mean(zap) #np.mean(np.mean(lindemann_per_atom(frames), axis=1))#lindemann_per_atom(frames)#   # type: ignore[no-any-return, no-untyped-call]
