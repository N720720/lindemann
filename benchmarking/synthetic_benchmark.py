import time

import numpy as np

from lindemann.index import per_atoms, per_frames, per_trj
from lindemann.index.mem_use import in_gb


def generate_test_data(num_frames, num_atoms):
    rng = np.random.default_rng(seed=42)
    return rng.random((num_frames, num_atoms, 3), dtype=np.float32)


def benchmark(function, positions, iterations=3):
    index = function(positions)  # Warm-up
    times = []
    for _ in range(iterations):
        start_time = time.time()
        function(positions)
        end_time = time.time()
        times.append(end_time - start_time)
    return np.mean(times), np.std(times), index


def main():
    num_frames = 5000
    num_atoms = 1103
    print(in_gb(num_frames, num_atoms))
    positions = generate_test_data(num_frames, num_atoms)
    iterations = 1

    mean_time, std_time, index = benchmark(per_trj.calculate, positions, iterations=iterations)
    print(f"Mean time of {iterations} runs: {mean_time:.6f} s ± {std_time:.6f} s")
    print(f"Index: {index}")


if __name__ == "__main__":
    main()
