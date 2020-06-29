import pytest

from lindemann.example import hello
from lindemann.index import per_atoms, per_frames, per_trj
from lindemann.trajectory import read
import numpy as np

"Testing the individal parts of the index module, its possible to change the test setup for individual modules"


@pytest.mark.parametrize(
    ("trajectory", "lindemannindex"),
    [
        ("tests/test_example/459_01.lammpstrj", np.float(0.025923892565654555)),
        ("tests/test_example/459_02.lammpstrj", np.float(0.026426709832984754)),
    ],
)
# def test_setup(trajectory: str, lindemannindex: np.float) -> bool:


def test_tra(trajectory, lindemannindex):
    """Example test with parametrization."""
    frame = read.frames(trajectory)
    assert per_trj.calculate(frame) == lindemannindex


@pytest.mark.parametrize(
    ("trajectory", "lindemannindex"),
    [
        ("tests/test_example/459_01.lammpstrj", np.float(0.025923892565654555)),
        ("tests/test_example/459_02.lammpstrj", np.float(0.026426709832984754)),
    ],
)
def test_frames(trajectory, lindemannindex):
    """Example test with parametrization."""
    frame = read.frames(trajectory)
    test_array = per_frames.calculate(frame)
    assert test_array[-1] == lindemannindex


@pytest.mark.parametrize(
    ("trajectory", "lindemannindex"),
    [
        ("tests/test_example/459_01.lammpstrj", np.float(0.025923892565654555)),
        ("tests/test_example/459_02.lammpstrj", np.float(0.026426709832984754)),
    ],
)
def test_atoms(trajectory, lindemannindex):
    """Example test with parametrization."""
    frame = read.frames(trajectory)
    assert np.mean(per_atoms.calculate(frame)[-1]) == lindemannindex
