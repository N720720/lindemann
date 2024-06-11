from typing import Any

import numba as nb
import numpy as np
import numpy.typing as npt


@nb.njit(fastmath=True, parallel=False)
def parallel_variance(
    n_a: npt.NDArray[np.float32],
    avg_a: npt.NDArray[np.float32],
    m2_a: npt.NDArray[np.float32],
    n_b: npt.NDArray[np.float32],
    avg_b: npt.NDArray[np.float32],
    m2_b: npt.NDArray[np.float32],
) -> tuple[npt.NDArray[np.float32], npt.NDArray[np.float32], npt.NDArray[np.float32]]:
    """
    Parallel variance computation using Welford's algorithm.

    Args:
        n_a (npt.NDArray[np.float32]): Count of samples in the first dataset.
        avg_a (npt.NDArray[np.float32]): Mean of the first dataset.
        m2_a (npt.NDArray[np.float32]): Second moment of the first dataset.
        n_b (npt.NDArray[np.float32]): Count of samples in the second dataset.
        avg_b (npt.NDArray[np.float32]): Mean of the second dataset.
        m2_b (npt.NDArray[np.float32]): Second moment of the second dataset.

    Returns:
        Tuple[npt.NDArray[np.float32], npt.NDArray[np.float32], npt.NDArray[np.float32]]: Combined count, mean, and second moment.
    """
    n_ab = n_a + n_b
    delta = avg_b - avg_a
    mean_ab = avg_a + delta * n_b / n_ab
    m2_ab = m2_a + m2_b + delta**2 * n_a * n_b / n_ab
    return n_ab, mean_ab, m2_ab


@nb.njit(fastmath=True, parallel=False)
def calculate_chunk(
    positions: npt.NDArray[np.float32], start_frame: int, end_frame: int
) -> tuple[npt.NDArray[np.float32], npt.NDArray[np.float32], int]:
    """
    Calculate the mean and variance for a chunk of frames.

    Args:
        positions (npt.NDArray[np.float32]): Array of shape (num_frames, num_atoms, 3) containing the positions.
        start_frame (int): The starting frame index for the chunk.
        end_frame (int): The ending frame index for the chunk.

    Returns:
        Tuple[npt.NDArray[np.float32], npt.NDArray[np.float32], int]: Mean distances, second moment distances, and count of frames processed.
    """
    num_frames, num_atoms, _ = positions.shape
    num_distances = num_atoms * (num_atoms - 1) // 2

    mean_distances = np.zeros(num_distances, dtype=np.float32)
    m2_distances = np.zeros(num_distances, dtype=np.float32)

    for frame in range(start_frame, end_frame):
        index = 0
        frame_count = frame + 1 - start_frame
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

    return mean_distances, m2_distances, frame_count


@nb.njit(fastmath=True, parallel=True)
def calculate(positions: npt.NDArray[np.float32], num_chunks: int) -> np.floating[Any]:
    """
    Calculate the Lindemann index in parallel using multiple chunks.

    Args:
        positions (npt.NDArray[np.float32]): Array of shape (num_frames, num_atoms, 3) containing the positions.
        num_chunks (int): Number of chunks to divide the frames into for parallel processing.

    Returns:
        float: The calculated Lindemann index.
    """
    num_frames, num_atoms, _ = positions.shape
    chunk_size = num_frames // num_chunks

    all_mean_distances = np.zeros((num_chunks, num_atoms * (num_atoms - 1) // 2), dtype=np.float32)
    all_m2_distances = np.zeros((num_chunks, num_atoms * (num_atoms - 1) // 2), dtype=np.float32)
    all_counts = np.zeros(num_chunks, dtype=np.int32)

    for chunk in nb.prange(num_chunks):
        start_frame = chunk * chunk_size
        end_frame = (chunk + 1) * chunk_size if chunk < num_chunks - 1 else num_frames
        mean_distances, m2_distances, count = calculate_chunk(positions, start_frame, end_frame)
        all_mean_distances[chunk] = mean_distances
        all_m2_distances[chunk] = m2_distances
        all_counts[chunk] = count

    final_mean_distances = all_mean_distances[0]
    final_m2_distances = all_m2_distances[0]
    final_count = all_counts[0]

    for chunk in range(1, num_chunks):
        final_count, final_mean_distances, final_m2_distances = parallel_variance(
            final_count,
            final_mean_distances,
            final_m2_distances,
            all_counts[chunk],
            all_mean_distances[chunk],
            all_m2_distances[chunk],
        )

    return np.mean(np.sqrt(final_m2_distances / final_count) / final_mean_distances)
