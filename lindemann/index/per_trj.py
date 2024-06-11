from typing import Any

import numba as nb
import numpy as np
import numpy.typing as npt

# @nb.njit(fastmath=True)
# def calculate(positions: npt.NDArray[np.float32]) -> np.floating[Any]:
#     """
#     Calculates the overall Lindemann index for a series of atomic positions over multiple frames.

#     Args:
#         positions (npt.NDArray[np.float32]): Array of atomic positions with shape (num_frames, num_atoms, 3).

#     Returns:
#         np.floating[np.Any]: The overall Lindemann index.
#     """
#     num_frames, num_atoms, _ = positions.shape
#     num_distances = num_atoms * (num_atoms - 1) // 2

#     mean_distances = np.zeros(num_distances, dtype=np.float32)
#     m2_distances = np.zeros(num_distances, dtype=np.float32)

#     for frame in range(num_frames):
#         index = 0
#         frame_count = frame + 1
#         for i in range(num_atoms):
#             for j in range(i + 1, num_atoms):
#                 dist = 0.0
#                 for k in range(3):
#                     dist += (positions[frame, i, k] - positions[frame, j, k]) ** 2
#                 dist = np.sqrt(dist)
#                 delta = dist - mean_distances[index]
#                 mean_distances[index] += delta / frame_count
#                 delta2 = dist - mean_distances[index]
#                 m2_distances[index] += delta * delta2

#                 index += 1

#     return np.mean(np.sqrt(m2_distances / num_frames) / mean_distances)


@nb.njit(fastmath=True)  # type: ignore
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
    # return np.mean(np.sum(np.nan_to_num(np.divide(np.sqrt(np.divide(array_var, nframes)), array_mean)), 1) / (natoms-1))
    return lindemann_indices


def calculate(frames: npt.NDArray[np.float64]) -> float:

    return np.mean(np.nanmean(lindemann_per_atom(frames), axis=1))  # type: ignore[no-any-return, no-untyped-call]
