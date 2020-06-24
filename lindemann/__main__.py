# type: ignore[attr-defined]

import random
from enum import Enum
from typing import Optional

import numpy as np
import typer
from rich.console import Console

from lindemann import __version__
from lindemann.example import hello
from lindemann.index.per_trj import lindemann_process_frames
from lindemann.trjread import frames


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
            f"[yellow]lindemann[/] version: [bold blue]{__version__}[/]"
        )
        raise typer.Exit()


@app.command(name="")
def main(
    trjfile: str = typer.Option(
        None, "-t", "--traj", help="The Trajectory file you want to process.",
    ),
    color: Optional[Color] = typer.Option(
        None,
        "-c",
        "--color",
        "--colour",
        case_sensitive=False,
        help="Color for trjfile. If not specified then choice will be random.",
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
    if color is None:
        # If no color specified use random value from `Color` class
        color = random.choice(list(Color.__members__.values()))
    trjfile: str = trjfile
    console.print(trjfile)
    my_frames = frames.frames(trjfile)
    indices = lindemann_process_frames(
        my_frames, len(my_frames), len(my_frames)
    )
    ind = np.mean(np.nanmean(indices, axis=1))
    console.print(ind)
    # index = lindemann.index.per_trj.lindemann_process_frames(trjfile)
    # greeting: float = index
    # console.print(f"The lindemann index is [bold {color}]{greeting}[/]")
