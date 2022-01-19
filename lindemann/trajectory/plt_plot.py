import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt

mpl.use("Agg")


def lindemann_vs_frames(indices: npt.NDArray[np.float64]) -> str:
    plt.figure(1)
    plt.title("Lindemann index per frame")
    plt.xlabel("Frames")
    plt.ylabel("Lindemann index")
    plt.plot(np.arange(0, len(indices)), indices, "+")
    plt.tight_layout()
    # plt.show()
    plt.savefig("lindemann_per_frame.pdf")
    return "lindemann_per_frame.pdf"
