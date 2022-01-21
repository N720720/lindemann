from typing import List

import numba as nb
import numpy as np
import numpy.typing as npt
from numba import float32


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


@nb.njit(fastmath=True, error_model="numpy")  # type: ignore # , cache=True) #(parallel=True)
def calculate(frames: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:

    """
    Calculate the lindemann index for each atom AND FRAME

    Return a ndarray of shape (len_frames, natoms, natoms)

    Warning this can produce extremly large ndarrays in memory 
    depending on the size of the cluster and the ammount of frames.
    """
    first = True
    # natoms = natoms
    dt = frames.dtype
    natoms = len(frames[0])
    nframes = len(frames)
    len_frames = len(frames)
    array_mean = np.zeros((natoms, natoms), dtype=dt)
    array_var = np.zeros((natoms, natoms), dtype=dt)
    # array_distance = np.zeros((natoms, natoms))
    iframe = dt.type(1)
    lindex_array = np.zeros((len_frames, natoms), dtype=dt)
    for q, coords in enumerate(frames):
        # print("processing frame {}/{}".format(iframe, nframes))
        # print(q)
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
        iframe += 1
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
            lindemann_indices = np.divide(np.sqrt(np.divide(array_var, nframes)), array_mean)
            lindemann_indices = np.asarray([np.mean(lin[lin != 0]) for lin in lindemann_indices])
            # print(type(lindemann_indices))
        # print(lindemann_indices.shape)
        lindex_array[q] = lindemann_indices  # np_mean(lindemann_indices,axis=1)
    return lindex_array
