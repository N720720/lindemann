import cProfile
import io
import pstats
from pathlib import Path

import numpy as np

from lindemann.main import main


def profile_lindemann(args, number_of_runs=5):

    main(args)  # Warmup
    all_times = []
    for run_number in range(number_of_runs):

        profiler = cProfile.Profile()

        profiler.enable()
        main(args)
        profiler.disable()

        stats_stream = io.StringIO()
        stats = pstats.Stats(profiler, stream=stats_stream).sort_stats("tottime")
        stats.print_stats()

        stats_stream.seek(0)

        for line in stats_stream:
            if "tottime" in line:
                line = next(stats_stream)
                time = float(line.split()[1])
                all_times.append(time)
                print(f"Time for run {run_number + 1}: {time} seconds")
                break

    all_times = np.array(all_times)
    mean_time = np.mean(all_times)
    std_time = np.std(all_times)
    print(f"Mean time of {number_of_runs} runs: {mean_time:.6f} s ± {std_time:.6f} s")


if __name__ == "__main__":
    args = [Path("tests/test_example/459_01.lammpstrj")]
    profile_lindemann(args, number_of_runs=5)
