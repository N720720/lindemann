import numpy as np
import numpy.typing as npt


def in_gb(frames: npt.NDArray[np.float64]) -> str:
    """Shows the size of the array in memory in GB.

    Args:
        frames (npt.NDArray[np.float64]): numpy array of shape(frames,atoms)

    Returns:
        str: Size of array in GB.
    """
    natoms = len(frames[0])
    nframes = len(frames)
    return f"This will use {np.round((np.zeros((natoms, natoms)).nbytes/1024**3),4)} GB"  # type: ignore[no-untyped-call]
