#!/usr/bin/env python

import os
import sys
import time

import numba as nb
import numpy as np
from ovito.io import import_file
from ovito.modifiers import SelectTypeModifier

# In[2]:

'''
@nb.njit(fastmath=True, cache=True)
def pdist_numba_py_sym(xs):
    """Unvectorized Python."""
    n, p = xs.shape
    A = np.zeros((n, n))
    for i in range(n):
        for j in range(i+1, n):
            d = 0.0
            for k in range(p):
                d += (xs[i, k] - xs[j, k])**2
            A[i,j] = np.sqrt(d)
    A += A.T
    return A
'''


@nb.njit(fastmath=True)  # , cache=True) #(parallel=True)
def lindemann_process_frames(frames, nframes, natoms):
    # natoms = natoms
    natoms = len(frames[0])
    array_mean = np.zeros((natoms, natoms))
    array_var = np.zeros((natoms, natoms))
    array_distance = np.zeros((natoms, natoms))
    iframe = 1
    for coords in frames:
        # print("processing frame {}/{}".format(iframe, nframes))

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
        if iframe > nframes:
            break

    for i in range(natoms):
        for j in range(i + 1, natoms):
            array_mean[j, i] = array_mean[i, j]
            array_var[j, i] = array_var[i, j]

    lindemann_indices = np.divide(
        np.sqrt(np.divide(array_var, nframes)), array_mean
    )
    return lindemann_indices


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

    natoms = len(frames)
    indices = lindemann_process_frames(frames, nframes, natoms)
    indices = np.mean(np.nanmean(indices, axis=1))
    print(indices)

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
