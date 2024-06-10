import numpy as np
import numpy.typing as npt


def in_gb(nframes: int, natoms: int) -> str:
    """
    Calculates and shows the size of the memory allocations related to
    the different flag options in gigabytes (GB).

    Args:
        nframes (int): The number of frames in the trajectory.
        natoms (int): The number of atoms per frame in the trajectory.

    Returns:
        str: A formatted string containing the memory usage for different configurations:
             - per_trj: Memory required when the `-t` flag is used.
             - per_frames: Memory required when the `-f` flag is used.
             - per_atoms: Memory required when the `-a` flag is used.

    This function assumes memory calculations based on numpy's float32 data type.
    """

    num_distances = natoms * (natoms - 1) // 2
    float_size = np.float32().nbytes
    trj = nframes * natoms * 3 * float_size
    atom_atom_array = 2 * natoms * natoms * float_size
    atom_array = natoms * float_size
    linde_index = nframes * natoms * float_size
    sum_bytes = trj + atom_atom_array + atom_array + linde_index
    per_trj = (
        f"\nFlag -t (per_trj) will use {np.round((trj+num_distances*2*float_size)/1024**3,4)} GB\n"
    )
    online_per_trj = (
        f"Flag -ot (per_trj) will use {np.round((num_distances*2*float_size)/1024**3,4)} GB\n"
    )
    per_frames = f"Flag -f (per_frames) will use {np.round((trj+(num_distances*2*float_size)+(nframes*float_size))/1024**3,4)} GB\n"
    online_per_frames = f"Flag -of (per_frames) will use {np.round(((num_distances*2*float_size)+(nframes*float_size))/1024**3,4)} GB\n"
    per_atoms = f"Flag -a (per_atoms) will use {np.round((sum_bytes)/1024**3,4)} GB\n"
    online_per_atoms = f"Flag -oa (per_atoms) will use {np.round((sum_bytes-trj)/1024**3,4)} GB"
    return f"{per_trj}{online_per_trj}{per_frames}{online_per_frames}{per_atoms}{online_per_atoms}"
