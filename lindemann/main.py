import re
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import typer
from psutil import cpu_count
from rich.console import Console

from lindemann import __version__
from lindemann.index import (
    mem_use,
    online_atoms,
    online_frames,
    online_trj,
    parallel_trj,
    per_atoms,
    per_frames,
    per_trj,
)
from lindemann.trajectory import plt_plot, read, save

app = typer.Typer(
    name="lindemann",
    help="""lindemann is a Python package to calculate the Lindemann index of a LAMMPS trajectory,
    as well as the progression of the Lindemann index per frame or per atom and frame of
    temperature ramps for phase transition analysis.""",
    add_completion=False,
)
console = Console()


def version_callback(value: bool):
    """Prints the version of the package."""
    if value:
        console.print(f"[magenta]lindemann[/] version: [bold blue]{__version__}[/]")
        raise typer.Exit()


@app.command()
def main(
    trjfile: list[Path] = typer.Argument(
        ...,
        help="The trajectory file(s). If no other option is selected, the Lindemann index is calculated for the trajectory. \
              Equivalent to the -t option. If you pass more than one trajectory they will be calculated in parallel. \
              Only works with no flag or -t flag.",
    ),
    trj: bool = typer.Option(
        False, "-t", help="Calculates the Lindemann-Index for the Trajectory file(s)"
    ),
    on_trj: bool = typer.Option(
        False,
        "-ot",
        help="Calculates the Lindemann-Index for the Trajectory file(s) (reduced memory usage)",
    ),
    par_trj: bool = typer.Option(
        False,
        "-pt",
        help="Calculates the Lindemann-Index for the Trajectory file(s) in parallel.",
    ),
    frames: bool = typer.Option(
        False, "-f", help="Calculates the Lindemann-Index for each frame."
    ),
    on_frames: bool = typer.Option(
        False, "-of", help="Calculates the Lindemann-Index for each frame. (reduced memory usage)"
    ),
    atoms: bool = typer.Option(
        False, "-a", help="Calculates the Lindemann-Index for each atom for each frame."
    ),
    on_atoms: bool = typer.Option(
        False,
        "-oa",
        help="Calculates the Lindemann-Index for each atom for each frame. (reduced memory usage)",
    ),
    plot: bool = typer.Option(False, "-p", help="Returns a plot Lindemann-Index vs. Frame."),
    lammpstrj: bool = typer.Option(
        False,
        "-l",
        help="Saves the individual Lindemann-Index of each Atom in a lammpstrj, so it can be viewed in Ovito.",
    ),
    version: bool = typer.Option(
        None,
        "-v",
        "--version",
        callback=version_callback,
        is_eager=True,
        help="Prints the version of the lindemann package.",
    ),
    timeit: bool = typer.Option(
        False, "-ti", "-timeit", help="Uses timeit module to show running time"
    ),
    mem_useage: bool = typer.Option(
        False,
        "-m",
        "-mem_use",
        help="Calculates the memory use. Run it before you use any of the CLI functionality despite the -t flag",
    ),
):
    """
    lindemann is a Python package to calculate the Lindemann index of a LAMMPS trajectory, as well
    as the progression of the Lindemann index per frame or per atom and frame of temperature ramps
    for phase transition analysis.
    """

    n_cores = min(len(trjfile), cpu_count())
    single_process = len(trjfile) == 1
    trjfile_str = [str(trjf) for trjf in trjfile]

    def calculate_single_pipeline(pipeline_func, data_func, save_filename=None, save_func=None):
        pipeline, data = pipeline_func(trjfile_str[0])
        results = data_func(pipeline, data)
        if save_filename and save_func:
            save_func(save_filename, results)
            console.print(f"[magenta]Lindemann index saved as:[/] [bold blue]{save_filename}[/]")
        else:
            console.print(
                f"[magenta]lindemann index for the Trajectory:[/] [bold blue]{results}[/]"
            )
        typer.Exit()

    def calculate_single(trjfile, calc_func, save_filename=None, save_func=None, cpu_count=None):
        frames = read.frames(trjfile)
        if cpu_count:
            results = calc_func(frames, cpu_count)
        else:
            results = calc_func(frames)
        if save_filename and save_func:
            save_func(save_filename, results)
            console.print(f"[magenta]Lindemann index saved as:[/] [bold blue]{save_filename}[/]")
        else:
            console.print(
                f"[magenta]lindemann index for the Trajectory:[/] [bold blue]{results}[/]"
            )
        typer.Exit()

    def calculate_parallel(trjfile_str, calc_func):
        trj_frames = [read.frames(tf) for tf in trjfile_str]
        with Pool(n_cores) as p:
            console.print(f"Using {n_cores} cores")
            res = p.map(calc_func, trj_frames)
            console.print(res)
        typer.Exit()

    if on_trj and single_process:
        calculate_single_pipeline(read.trajectory, online_trj.calculate)
    elif trj and single_process:
        calculate_single(trjfile_str[0], per_trj.calculate)
    elif par_trj and single_process:
        calculate_single(trjfile_str[0], parallel_trj.calculate, cpu_count=cpu_count())
    elif trj and not single_process:
        calculate_parallel(trjfile_str, per_trj.calculate)
    elif frames and single_process:
        calculate_single(
            trjfile_str[0], per_frames.calculate, "lindemann_index_per_frame.txt", np.savetxt
        )
    elif frames and not single_process:
        console.print("multiprocessing is implemented only for the -t flag")
        typer.Exit()
    elif on_frames and single_process:
        calculate_single_pipeline(
            read.trajectory, online_frames.calculate, "lindemann_index_per_frame.txt", np.savetxt
        )
    elif on_frames and not single_process:
        console.print("multiprocessing is implemented only for the -t flag")
        typer.Exit()
    elif atoms and single_process:
        calculate_single(
            trjfile_str[0], per_atoms.calculate, "lindemann_index_per_atom.txt", np.savetxt
        )
    elif atoms and not single_process:
        console.print("multiprocessing is implemented only for the -t flag")
        typer.Exit()
    elif on_atoms and single_process:
        calculate_single_pipeline(
            read.trajectory, online_atoms.calculate, "lindemann_index_per_atoms.txt", np.savetxt
        )
    elif on_atoms and not single_process:
        console.print("multiprocessing is implemented only for the -t flag")
        typer.Exit()
    elif plot and single_process:
        tjr_frames = read.frames(trjfile_str[0])
        indices = per_frames.calculate(tjr_frames)
        plot_filename = plt_plot.lindemann_vs_frames(indices)
        console.print(f"[magenta]Saved file as:[/] [bold blue]{plot_filename}[/]")
        typer.Exit()
    elif plot and not single_process:
        console.print("multiprocessing is implemented only for the -t flag")
        typer.Exit()
    elif lammpstrj and single_process:
        calculate_single(trjfile_str[0], per_atoms.calculate, trjfile_str[0], save.to_lammps)
        typer.Exit()
    elif lammpstrj and not single_process:
        console.print("multiprocessing is implemented only for the -t flag")
        typer.Exit()
    elif timeit and single_process:
        tjr_frames = read.frames(trjfile_str[0])
        start = time.time()
        linde_for_time = per_trj.calculate(tjr_frames)
        time_diff = time.time() - start
        console.print(
            f"[magenta]lindemann index for the Trajectory:[/] [bold blue]{linde_for_time}[/] \n"
            f"[magenta]Runtime:[/] [bold green]{time_diff}[/]"
        )
        typer.Exit()
    elif timeit and not single_process:
        console.print("multiprocessing is implemented only for the -t flag")
        typer.Exit()
    elif mem_useage and single_process:
        tjr_frames = read.frames(trjfile_str[0])
        nframes, natoms, _ = tjr_frames.shape
        mem_use_in_gb = mem_use.in_gb(nframes, natoms)
        console.print(f"[magenta]Memory use:[/] [bold blue]{mem_use_in_gb}[/]")
        typer.Exit()
    elif mem_useage and not single_process:
        console.print("multiprocessing is implemented only for the -t flag")
        typer.Exit()
    else:
        if single_process:
            calculate_single_pipeline(read.trajectory, online_trj.calculate)
        else:
            calculate_parallel(trjfile_str, per_trj.calculate)


if __name__ == "__main__":
    app()
