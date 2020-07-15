import numpy as np


def in_gb(frames: np.ndarray) -> str:
    natoms = len(frames[0])
    nframes = len(frames)
    return f"This will use {np.round((np.zeros((nframes, natoms, natoms)).nbytes/1024**3),4)} GB"
