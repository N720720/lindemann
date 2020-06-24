import os

import numpy as np
from ovito.io import import_file
from ovito.modifiers import SelectTypeModifier


def frames(trjfile, nframes=None):

    """
    Get the frames from the lammps trajectory using ovito pipeline and import_file function.
    It returns frames and the number of frames to use for calculating the Lindemann Index.
    """

    if not os.path.exists(trjfile):
        raise RuntimeError(f"Error: file {trjfile} not found!")

    pipeline = import_file(trjfile, sort_particles=True)
    num_frame = pipeline.source.num_frames
    pipeline.modifiers.append(
        SelectTypeModifier(
            operate_on="particles", property="Particle Type", types={1, 2, 3}
        )
    )
    data = pipeline.compute()
    num_particle = data.particles.count

    # If no argument is given use all frames
    if nframes is None:
        nframes = num_frame
    # make sure nobody puts more frames then exists
    assert num_frame >= nframes

    # initialise array, could be problematic for big clusters and a lot of frames
    position = np.zeros((nframes, num_particle, 3))

    for frame in range(nframes):
        data = pipeline.compute(frame)
        position[frame, :, :] = np.array(data.particles["Position"])
    frames = position

    return frames
