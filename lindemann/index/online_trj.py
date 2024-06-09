from typing import Optional

import numba as nb
import numpy as np
from ovito.data import DataCollection
from ovito.pipeline import Pipeline

from lindemann.trajectory import read


@nb.njit(fastmath=True, parallel=False)
def calculate_frame(positions, mean_distances, m2_distances, frame: int, num_atoms: int):

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


def calculate(pipeline: Pipeline, data: DataCollection, nframes: Optional[int] = None):
    num_particle = data.particles.count
    num_frame = pipeline.source.num_frames
    if nframes is None:
        nframes = num_frame
    elif nframes > num_frame:
        raise ValueError(f"Requested {nframes} frames, but only {num_frame} frames are available.")

    num_distances = num_particle * (num_particle - 1) // 2
    mean_distances = np.zeros(num_distances, dtype=np.float32)
    m2_distances = np.zeros(num_distances, dtype=np.float32)
    print(nframes, num_particle)
    for frame in range(nframes):
        data = pipeline.compute(frame)
        calculate_frame(
            data.particles["Position"].array, mean_distances, m2_distances, frame, num_particle
        )

    linde = np.mean(np.sqrt(m2_distances / nframes) / mean_distances)
    return linde
