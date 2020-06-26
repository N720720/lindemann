import numba as nb
import numpy as np


@nb.njit(fastmath=True)  # , cache=True) #(parallel=True)
def lindemann_per_frames(frames: np.ndarray) -> np.ndarray:

    """
    Calculate the lindemann index for each atom AND FRAME

    Return a ndarray of shape (len_frames, natoms, natoms)

    Warning this can produce extremly large ndarrays in memory 
    depending on the size of the cluster and the ammount of frames.
    """
    # natoms = natoms
    natoms = len(frames[0])
    nframes = len(frames)
    len_frames = len(frames)
    array_mean = np.zeros((natoms, natoms))
    array_var = np.zeros((natoms, natoms))
    # array_distance = np.zeros((natoms, natoms))
    iframe = 1
    lindex_array = np.zeros((len_frames, natoms, natoms))
    for q, coords in enumerate(frames):
        # print("processing frame {}/{}".format(iframe, nframes))
        # print(q)
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

        lindemann_indices = np.divide(
            np.sqrt(np.divide(array_var, nframes)), array_mean
        )
        # lindemann_indices = np.nanmean(np.sqrt(array_var/nframes)/array_mean, axis=1)
        lindex_array[q] = lindemann_indices
    return lindex_array


def calculate(indices: np.ndarray) -> np.ndarray:
    """
    Small helper function, since numba has not implemented the np.nanmean with axis parameter 
    I cant implemnet this in the jit function for now.
    """
    return [np.nanmean(i, axis=1) for i in lindemann_per_frames(indices)]
