#!/usr/bin/env python

import os
import sys
import time

import numba as nb
import numpy as np
from ovito.io import import_file
from ovito.modifiers import SelectTypeModifier


def process_trjfile(trjfile, nframes=None):

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
    position = np.zeros((num_frame, num_particle, 3))
    for frame in range(num_frame):
        data = pipeline.compute(frame)
        position[frame, :, :] = np.array(data.particles["Position"])
    frames = position
    # print(frames.shape)
    if nframes is None:
        nframes = len(frames)
    assert len(frames) >= nframes

    return frames, nframes
