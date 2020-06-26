# type: ignore[attr-defined]

import random
from enum import Enum
from typing import Optional

import numpy as np
import typer
from rich.console import Console

from lindemann import __version__
from lindemann.example import hello

# from lindemann.index.per_trj import lindemann_process_frames
# from lindemann.trjread import frames

app = typer.Typer(
    name="lindemann",
    help="lindemann is a python package to calculate the Lindemann index  of a lammps trajectory as well as the progression of the Lindemann index per frame of temperature ramps  for phase transition analysis.",
    add_completion=False,
)
console = Console()


def version_callback(value: bool):
    """Prints the version of the package."""
    if value:
        console.print(
            f"[magenta]lindemann[/] version: [bold blue]{__version__}[/]"
        )
        raise typer.Exit()


@app.command(name="")
def main(
    trjfile: str = typer.Argument(
        ...,
        # help="The Lammps Trajectory File, you can give the filepath or the filename if you are in the directory",
    ),
    trj: bool = typer.Option(
        False, help="Calculates the Lindemann-Index for the Trajectory file"
    ),
    frames: bool = typer.Option(
        False, help="Calculates the Lindemann-Index for each frame."
    ),
    atoms: bool = typer.Option(
        False,
        help="Calculates the Lindemann-Index for each atom for each frame.",
    ),
    plot: bool = typer.Option(
        False, help="Returns a plot Lindemann-Index vs. Frame."
    ),
    lammpstrj: bool = typer.Option(
        False,
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
):

    """
    lindemann is a python package to calculate the Lindemann index  of a lammps trajectory as well as the progression of the Lindemann index per frame of temperature ramps  for phase transition analysis.
    """


'''
class Color(str, Enum):
    white = "white"
    red = "red"
    cyan = "cyan"
    magenta = "magenta"
    yellow = "yellow"
    green = "green"


app = typer.Typer(
    name="lindemann",
    help="lindemann is a python package to calculate the Lindemann index  of a lammps trajectory as well as the progression of the Lindemann index per frame of temperature ramps  for phase transition analysis.",
    add_completion=False,
)
console = Console()


def version_callback(value: bool):
    """Prints the version of the package."""
    if value:
        console.print(
            f"[magenta]lindemann[/] version: [bold blue]{__version__}[/]"
        )
        raise typer.Exit()


@app.command(name="")
def main(
    trjfile: str = typer.Argument(
        None, help="The Trajectory file you want to process.",
    ),
    version: bool = typer.Option(
        None,
        "-v",
        "--version",
        callback=version_callback,
        is_eager=True,
        help="Prints the version of the lindemann package.",
    ),
    index: bool = typer.Option(
        None,
        "-i",
        "--index",
        help="Calculates the lindemann index for the given lammps trajectory.",
    ),
):
    """Prints a greeting for a giving trjfile."""
    trjfile: str = trjfile

    console.print(trjfile)
    # my_frames = frames.frames(trjfile)
    # indices = lindemann_process_frames(
    #    my_frames, len(my_frames), len(my_frames)
    # )
    # ind = np.mean(np.nanmean(indices, axis=1))
    # console.print(ind)
    # index = lindemann.index.per_trj.lindemann_process_frames(trjfile)
    # greeting: float = index
    # console.print(f"The lindemann index is [bold {color}]{greeting}[/]")
'''
