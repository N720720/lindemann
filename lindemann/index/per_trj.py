import numba as nb
import numpy as np


@nb.njit(fastmath=True)
def calculate(positions):
    num_frames, num_atoms, _ = positions.shape
    num_distances = num_atoms * (num_atoms - 1) // 2

    mean_distances = np.zeros(num_distances, dtype=np.float32)
    m2_distances = np.zeros(num_distances, dtype=np.float32)

    for frame in range(num_frames):
        index = 0
        frame_count = frame + 1
        for i in range(num_atoms):
            for j in range(i + 1, num_atoms):
                dist = 0.0
                for k in range(3):
                    dist += (positions[frame, i, k] - positions[frame, j, k]) ** 2
                dist = np.sqrt(dist)
                delta = dist - mean_distances[index]
                mean_distances[index] += delta / frame_count
                delta2 = dist - mean_distances[index]
                m2_distances[index] += delta * delta2

                index += 1

    return np.mean(np.sqrt(m2_distances / num_frames) / mean_distances)
