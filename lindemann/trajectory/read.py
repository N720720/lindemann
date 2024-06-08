from typing import Optional

import numpy as np
import numpy.typing as npt
from ovito.io import import_file
from ovito.modifiers import SelectTypeModifier


def frames(trjfile: str, nframes: Optional[int] = None) -> npt.NDArray[np.float32]:
    """
    Extracts the frame position data from a MD trajectory file using the OVITO pipeline.

    The function loads the specified trajectory file, applies a selection modifier to filter
    particles of type 1, 2, and 3, and computes the positions for a specified number of frames.
    If `nframes` is None, the function will attempt to process all frames in the trajectory.

    Parameters:
        trjfile (str): Path to the trajectory file to be processed.
        nframes (Optional[int]): The number of frames to process. If not specified, all frames
                                 in the trajectory file are processed. If the specified number
                                 exceeds the available frames in the file, a ValueError is raised.

    Returns:
        npt.NDArray[np.float32]: A 3D NumPy array of shape (nframes, num_particles, 3) containing
                                 the position data for each particle across the specified frames.

    Raises:
        ValueError: If `nframes` is more than the number of available frames in the trajectory file.

    Example:
        >>> positions = frames("path/to/trajectory.lammpstrj", 100)
        This would load 100 frames from the specified file and return the position data.
    """

    pipeline = import_file(trjfile, sort_particles=True)
    num_frame = pipeline.source.num_frames
    pipeline.modifiers.append(
        SelectTypeModifier(operate_on="particles", property="Particle Type", types={1, 2, 3})
    )
    data = pipeline.compute()
    num_particle = data.particles.count

    if nframes is None:
        nframes = num_frame
    elif nframes > num_frame:
        raise ValueError(f"Requested {nframes} frames, but only {num_frame} frames are available.")

    position = np.zeros((nframes, num_particle, 3), dtype=np.float32)

    for frame in range(nframes):
        data = pipeline.compute(frame)
        position[frame, :, :] = np.array(data.particles["Position"])
    frames = position

    return frames
