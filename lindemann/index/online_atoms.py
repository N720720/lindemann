from typing import Optional

import numba as nb
import numpy as np
import numpy.typing as npt
from ovito.data import DataCollection
from ovito.pipeline import Pipeline


@nb.njit(fastmath=True, error_model="numpy")  # type: ignore # , cache=True) #(parallel=True)
def calculate_frame(positions, array_mean, array_var, frame, natoms) -> npt.NDArray[np.float32]:
    """Calculate the contribution of each atom to the lindemann index over the frames

    Args:
        frames: numpy array of shape(frames,atoms)
    Returns:
        npt.NDArray[np.float32]: Returns 1D array with the progression of the lindeman index per frame of shape(frames, atoms)
    """

    frame_count = frame + 1
    for i in range(natoms):
        for j in range(i + 1, natoms):
            dist = 0.0
            for k in range(3):
                dist += (positions[i, k] - positions[j, k]) ** 2
            dist = np.sqrt(dist)
            mean = array_mean[i, j]
            var = array_var[i, j]
            delta = dist - mean
            update_mean = mean + delta / frame_count
            array_mean[i, j] = update_mean
            array_mean[j, i] = update_mean
            delta2 = dist - array_mean[i, j]
            update_var = var + delta * delta2
            array_var[i, j] = update_var
            array_var[j, i] = update_var

    np.fill_diagonal(array_mean, 1.0)
    lindemann_indices = np.divide(np.sqrt(np.divide(array_var, frame_count)), array_mean).astype(
        np.float32
    )
    lindemann_indices = np.asarray(
        [np.nanmean(lin[lin != 0]) for lin in lindemann_indices]
    ).astype(np.float32)

    return lindemann_indices


def calculate(pipeline: Pipeline, data: DataCollection, nframes: Optional[int] = None):
    num_particle = data.particles.count
    num_frame = pipeline.source.num_frames
    if nframes is None:
        nframes = num_frame
    elif nframes > num_frame:
        raise ValueError(f"Requested {nframes} frames, but only {num_frame} frames are available.")

    array_mean = np.zeros((num_particle, num_particle), dtype=np.float32)
    array_var = np.zeros((num_particle, num_particle), dtype=np.float32)
    lindex_array = np.zeros((nframes, num_particle), dtype=np.float32)
    for frame in range(nframes):
        data = pipeline.compute(frame)
        lindex_array[frame] = calculate_frame(
            data.particles["Position"].array, array_mean, array_var, frame, num_particle
        )
    return lindex_array
