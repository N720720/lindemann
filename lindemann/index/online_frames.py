from typing import Any, Optional

import numba as nb
import numpy as np
import numpy.typing as npt
from ovito.data import DataCollection
from ovito.pipeline import Pipeline


@nb.njit(fastmath=True, parallel=False)
def calculate_frame(
    positions: npt.NDArray[np.float32],
    mean_distances: npt.NDArray[np.float32],
    m2_distances: npt.NDArray[np.float32],
    frame: int,
    num_atoms: int,
) -> np.floating[Any]:
    """
    Calculates the Lindemann Index for a specific frame.

    Args:
        positions (npt.NDArray[np.float32]): Array of atomic positions for the current frame.
        mean_distances (npt.NDArray[np.float32]): Array to store the mean distances between pairs of atoms.
        m2_distances (npt.NDArray[np.float32]): Array to store the squared differences of distances between pairs of atoms.
        frame (int): The current frame index.
        num_atoms (int): The number of atoms.

    Returns:
        float: The Lindemann index for the current frame.
    """
    index = 0
    frame_count = frame + 1
    for i in range(num_atoms):
        for j in range(i + 1, num_atoms):
            dist = 0.0
            for k in range(3):
                dist += (positions[i, k] - positions[j, k]) ** 2

            dist = np.sqrt(dist)
            delta = dist - mean_distances[index]
            mean_distances[index] += delta / frame_count
            delta2 = dist - mean_distances[index]
            m2_distances[index] += delta * delta2

            index += 1
    return np.mean(np.sqrt(m2_distances / frame_count) / mean_distances)


def calculate(
    pipeline: Pipeline, data: DataCollection, nframes: Optional[int] = None
) -> npt.NDArray[np.float32]:
    """
    Calculates the Lindemann indices for a series of frames from an OVITO pipeline.

    Args:
        pipeline (Pipeline): The OVITO pipeline object.
        data (DataCollection): The data collection object from OVITO.
        nframes (Optional[int]): The number of frames to process. If None, all frames are processed.

    Returns:
        npt.NDArray[np.float32]: Array of Lindemann indices for each frame.

    Raises:
        ValueError: If the requested number of frames exceeds the available frames in the pipeline.
    """
    num_particle = data.particles.count
    num_frame = pipeline.source.num_frames
    if nframes is None:
        nframes = num_frame
    elif nframes > num_frame:
        raise ValueError(f"Requested {nframes} frames, but only {num_frame} frames are available.")

    num_distances = num_particle * (num_particle - 1) // 2
    mean_distances = np.zeros(num_distances, dtype=np.float32)
    m2_distances = np.zeros(num_distances, dtype=np.float32)
    lindemann_index_array = np.zeros(nframes, dtype=np.float32)
    for frame in range(nframes):
        data = pipeline.compute(frame)
        lindemann_index_array[frame] = calculate_frame(
            data.particles["Position"].array, mean_distances, m2_distances, frame, num_particle
        )
    return lindemann_index_array
