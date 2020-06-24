#!/usr/bin/env python

# In[1]:


import os
import sys
import time

import numba as nb
import numpy as np
import ovito
from ovito.io import export_file, import_file
from ovito.modifiers import SelectTypeModifier

# In[2]:


@nb.njit(fastmath=True)  # , cache=True) #(parallel=True)
def lindemann_process_frames(frames, nframes, natoms):
    # natoms = natoms
    natoms = len(frames[0])
    len_frames = len(frames)
    array_mean = np.zeros((natoms, natoms))
    array_var = np.zeros((natoms, natoms))
    # array_distance = np.zeros((natoms, natoms))
    iframe = 1
    lindex_array = np.zeros((len_frames, natoms, natoms))
    for q, coords in enumerate(frames):
        # print("processing frame {}/{}".format(iframe, nframes))
        print(q)
        n, p = coords.shape
        array_distance = np.zeros((n, n))
        for i in range(n):
            for j in range(i + 1, n):
                d = 0.0
                for k in range(p):
                    d += (coords[i, k] - coords[j, k]) ** 2
                array_distance[i, j] = np.sqrt(d)
        array_distance += array_distance.T

        #################################################################################
        # update mean and var arrays based on Welford algorithm suggested by Donald Knuth
        #################################################################################
        for i in range(natoms):
            for j in range(i + 1, natoms):
                xn = array_distance[i, j]
                mean = array_mean[i, j]
                var = array_var[i, j]
                delta = xn - mean
                # update mean
                array_mean[i, j] = mean + delta / iframe
                # update variance
                array_var[i, j] = var + delta * (xn - array_mean[i, j])
        iframe += 1
        if iframe > nframes + 1:
            break

        for i in range(natoms):
            for j in range(i + 1, natoms):
                array_mean[j, i] = array_mean[i, j]
                array_var[j, i] = array_var[i, j]

        lindemann_indices = np.divide(
            np.sqrt(np.divide(array_var, nframes)), array_mean
        )
        # lindemann_indices = np.nanmean(np.sqrt(array_var/nframes)/array_mean, axis=1)
        lindex_array[q] = lindemann_indices
    return lindex_array


# In[6]:


def process_trjfile(trjfile, nframes=None):
    """parse atom coordinates for all trajectories"""

    if not os.path.exists(trjfile):
        raise RuntimeError(f"Error: file {trjfile} not found!")

    #####################################################################################
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
    ####################################################################################

    # print(len(natoms))
    print(f"frames: {len(frames)}")
    print(nframes)
    natoms = len(frames[0])
    print(f"natoms: {natoms}")
    print(
        f"this will use {np.round((np.zeros((len(frames), natoms, natoms)).nbytes/1024**3),4)} GB"
    )
    indices = lindemann_process_frames(frames, nframes, natoms)
    print(indices.shape)
    print(indices[1])
    indices_per_atom = [np.nanmean(i, axis=1) for i in indices]
    indices = [np.mean(np.nanmean(i, axis=1)) for i in indices]
    print(indices_per_atom[-1])
    print(indices)

    pipeline = import_file(trjfile, sort_particles=True)
    # liste = []
    for frame, linde in enumerate(indices_per_atom):
        data = pipeline.compute(frame)
        data.particles_.create_property("lindemann", data=linde)
        export_file(
            data,
            f"outputfile{frame}.dump",
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
    # export_file(ldata, "output_data.txt",
    #            format = "lammps/dump",
    #            columns = ["Frame","Position.X", "Position.Y", "Position.Z","lindemann"],
    #            multiple_frames = True
    #            )
    # data = pipeline.compute()
    # pipeline = import_file(trjfile, sort_particles=True)
    # data = data.particles.create_property('lindemann', data=indices_per_atom[0])
    # export_file(data, "outputfile.dump", "lammps/dump",columns = ["Particle Identifier", "Particle Type", "Position.X", "Position.Y", "Position.Z","lindemann"],multiple_frames = True)
    filenames = [
        f"outputfile{frame}.dump" for frame in range(len(indices_per_atom))
    ]
    with open("myfile.lammpstrj", "w") as outfile:
        for fname in filenames:
            with open(fname) as infile:
                outfile.write(infile.read())
    for frame in range(len(indices_per_atom)):
        os.remove(f"outputfile{frame}.dump")
    return frames


def main():
    import argparse

    version = (
        "%(prog)s " + "<<VERSION>>" + "; updated at: " + "<<UPDATED-AT()>>"
    )
    desc = "default description"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("-v", "--version", version=version, action="version")
    parser.add_argument("trjfile", help="lammpstrj filename")
    parser.add_argument(
        "-n",
        dest="nframes",
        type=int,
        help="the maxium frames will be processed",
    )

    if len(sys.argv) == 1:
        parser.print_help()
        return

    cmdl = parser.parse_args()
    process_trjfile(cmdl.trjfile, nframes=cmdl.nframes)


if __name__ == "__main__":
    start = time.time()
    main()
    time_diff = time.time() - start
    print(time_diff)
# In[ ]:


# [0.01715329 0.0161475  0.0163233 ... 0.01074006 0.0164681  0.01535047]


# In[ ]:
