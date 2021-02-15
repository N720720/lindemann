import numba as nb
import numpy as np


@nb.njit(fastmath=True)  # , cache=True)
# def lindemann_per_atom(frames: np.ndarray) -> np.ndarray:
def lindemann_per_atom(frames):

    """Calculates the lindemann index for """
    natoms = len(frames[0])
    nframes = len(frames)
    array_mean = np.zeros((natoms, natoms))
    array_var = np.zeros((natoms, natoms))
    array_distance = np.zeros((natoms, natoms))
    iframe = 1
    for coords in frames:

        # here we do someting similar to scipy's spatial.distance.pdist scipy.spatial.distance.pdist unfortunaltelly
        n, p = coords.shape
        array_distance = np.zeros((n, n))
        for i in range(n):
            for j in range(i + 1, n):
                d = 0.0
                for k in range(p):
                    d += (coords[i, k] - coords[j, k]) ** 2
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
        iframe += 1
        if iframe > nframes:
            break

    for i in range(natoms):
        for j in range(i + 1, natoms):
            array_mean[j, i] = array_mean[i, j]
            array_var[j, i] = array_var[i, j]

    lindemann_indices = np.divide(
        np.sqrt(np.divide(array_var, nframes)), array_mean
    )
    return lindemann_indices


# def calculate(frames: np.ndarray) -> float:
def calculate(frames):
    """
    Small helper function, since numba has not implemented the np.nanmean with axis parameter 
    I cant implemnet this in the jit function for now.
    """

    indices = lindemann_per_atom(frames)
    index = np.mean(np.nanmean(indices, axis=1))

    return index
