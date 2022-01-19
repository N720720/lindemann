import numpy as np
import pytest

from lindemann.index import per_atoms, per_frames, per_trj
from lindemann.trajectory import read

"Testing the individal parts of the index module, its possible to change the test setup for individual modules"


@pytest.mark.parametrize(
    ("trajectory", "lindemannindex"),
    [
        ("tests/test_example/459_01.lammpstrj", np.round(float(0.025923892565654555), 12),),
        ("tests/test_example/459_02.lammpstrj", np.round(float(0.026426709832984754), 12),),
    ],
)
# def test_setup(trajectory: str, lindemannindex: float) -> bool:


def test_tra(trajectory, lindemannindex):
    """Example test with parametrization."""
    frame = read.frames(trajectory)
    assert np.round(per_trj.calculate(frame), 12) == lindemannindex


@pytest.mark.parametrize(
    ("trajectory", "lindemannindex"),
    [
        ("tests/test_example/459_01.lammpstrj", np.round(float(0.025923892565654555), 12),),
        ("tests/test_example/459_02.lammpstrj", np.round(float(0.026426709832984754), 12),),
    ],
)
def test_frames(trajectory, lindemannindex):
    """Example test with parametrization."""
    frame = read.frames(trajectory)
    test_array = per_frames.calculate(frame)
    assert np.round(test_array[-1], 12) == lindemannindex


@pytest.mark.parametrize(
    ("trajectory", "lindemannindex"),
    [
        ("tests/test_example/459_01.lammpstrj", np.round(float(0.025923892565654555), 12),),
        ("tests/test_example/459_02.lammpstrj", np.round(float(0.026426709832984754), 12),),
    ],
)
def test_atoms(trajectory, lindemannindex):
    """Example test with parametrization."""
    frame = read.frames(trajectory)
    assert np.round(np.mean(per_atoms.calculate(frame)[-1]), 12) == lindemannindex
