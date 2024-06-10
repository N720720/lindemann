import numba as nb
import numpy as np
import numpy.typing as npt


@nb.njit(fastmath=True, error_model="numpy")  # type: ignore # , cache=True) #(parallel=True)
def calculate(frames: npt.NDArray[np.float32]) -> npt.NDArray[np.float32]:
    """Calculate the contribution of each atom to the lindemann index over the frames

    Args:
        frames: numpy array of shape(frames,atoms)
    Returns:
        npt.NDArray[np.float32]: Returns 1D array with the progression of the lindeman index per frame of shape(frames, atoms)
    """
    len_frames, natoms, _ = frames.shape

    array_mean = np.zeros((natoms, natoms), dtype=np.float32)
    array_var = np.zeros((natoms, natoms), dtype=np.float32)
    lindex_array = np.zeros((len_frames, natoms), dtype=np.float32)
    for frame, coords in enumerate(frames):
        frame_count = frame + 1
        for i in range(natoms):
            for j in range(i + 1, natoms):
                dist = 0.0
                for k in range(3):
                    dist += (coords[i, k] - coords[j, k]) ** 2
                dist = np.sqrt(dist)
                mean = array_mean[i, j]
                var = array_var[i, j]
                delta = dist - mean
                update_mean = mean + delta / frame_count
                array_mean[i, j] = update_mean
                array_mean[j, i] = update_mean
                delta2 = dist - array_mean[i, j]
                update_var = var + delta * delta2
                array_var[i, j] = update_var
                array_var[j, i] = update_var

        np.fill_diagonal(array_mean, 1.0)
        lindemann_indices = np.divide(
            np.sqrt(np.divide(array_var, frame_count)), array_mean
        ).astype(np.float32)
        lindemann_indices = np.asarray(
            [np.nanmean(lin[lin != 0]) for lin in lindemann_indices]
        ).astype(np.float32)

        lindex_array[frame] = lindemann_indices
    return lindex_array
