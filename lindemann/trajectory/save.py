import os

import numpy as np
import numpy.typing as npt
from ovito.io import export_file, import_file
from ovito.modifiers import SelectTypeModifier

"""
I had big problems to get the ovito export_file module to do what I want, to put a numpy ndarry to a lammpstrj.
"""


def to_lammps(trjfile: str, indices_per_atom: npt.NDArray[np.float64]) -> str:
    pipeline = import_file(trjfile, sort_particles=True)

    for frame, linde in enumerate(indices_per_atom):
        data = pipeline.compute(frame)
        data.particles_.create_property("lindemann", data=linde)
        export_file(
            data,
            f"lindemann_outputfile_X{frame}.dump",
            "lammps/dump",
            columns=[
                "Particle Identifier",
                "Particle Type",
                "Position.X",
                "Position.Y",
                "Position.Z",
                "lindemann",
            ],
            multiple_frames=True,
        )

    """
    This is section is weird and very hacky, I dont like it. Does someone have a idea how to do this in a better way?
    First I save the files with the ovito export_file, then I save them in a file, finally I remove the files... I know...
    """

    filenames = [f"lindemann_outputfile_X{frame}.dump" for frame in range(len(indices_per_atom))]

    file_name = "lindemann_per_atom.lammpstrj"
    with open(file_name, "w") as outfile:
        for fname in filenames:
            with open(fname) as infile:
                outfile.write(infile.read())
    for frame in range(len(indices_per_atom)):
        os.remove(f"lindemann_outputfile_X{frame}.dump")
    return "saved trajectory as lindemann_per_atom.lammpstrj"
