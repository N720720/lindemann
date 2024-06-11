from typing import Optional

import numba as nb
import numpy as np
import numpy.typing as npt
from ovito.data import DataCollection
from ovito.pipeline import Pipeline


@nb.njit(fastmath=True, parallel=False)
def calculate_frame(
    positions: npt.NDArray[np.float32],
    array_mean: npt.NDArray[np.float32],
    array_var: npt.NDArray[np.float32],
    frame: int,
    natoms: int,
) -> npt.NDArray[np.float32]:
    """
    Calculates the contribution of the individual atomic positions to the Lindemann Index for a specific frame.

    Args:
        positions (npt.NDArray[np.float32]): Array of atomic positions for the current frame.
        array_mean (npt.NDArray[np.float32]): Array to store the mean distances.
        array_var (npt.NDArray[np.float32]): Array to store the variance of distances.
        frame (int): The current frame index.
        natoms (int): The number of atoms.

    Returns:
        npt.NDArray[np.float32]: Array of the individual atomic contributions to the Lindemann indices for the current frame.
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


def calculate(
    pipeline: Pipeline, data: DataCollection, nframes: Optional[int] = None
) -> npt.NDArray[np.float32]:
    """
    Calculates the contribution of the individual atomic positions to the Lindemann Index for a series of frames from an OVITO pipeline.

    Args:
        pipeline (Pipeline): The OVITO pipeline object.
        data (DataCollection): The data collection object from OVITO.
        nframes (Optional[int]): The number of frames to process. If None, all frames are processed.

    Returns:
        npt.NDArray[np.float32]: Array of the individual atomic contributions to the Lindemann indices for each frame.

    Raises:
        ValueError: If the requested number of frames exceeds the available frames in the pipeline.
    """
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
