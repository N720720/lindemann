import numpy as np
import pytest

from lindemann.index import online_frames, online_trj, per_atoms, per_frames, per_trj
from lindemann.trajectory import read

"Testing the individal parts of the index module, its possible to change the test setup for individual modules"


@pytest.mark.parametrize(
    ("trajectory", "lindemannindex"),
    [
        (
            "tests/test_example/459_01.lammpstrj",
            0.025923892565654555,
        ),
        (
            "tests/test_example/459_02.lammpstrj",
            0.026426709832984754,
        ),
    ],
)
def test_tra(trajectory, lindemannindex):
    """Example test with parametrization."""
    frame = read.frames(trajectory)
    assert np.isclose(per_trj.calculate(frame), lindemannindex)


@pytest.mark.parametrize(
    ("trajectory", "lindemannindex"),
    [
        (
            "tests/test_example/459_01.lammpstrj",
            0.025923892565654555,
        ),
        (
            "tests/test_example/459_02.lammpstrj",
            0.026426709832984754,
        ),
    ],
)
def test_online_tra(trajectory, lindemannindex):
    """Example test with parametrization."""
    pipeline, data = read.trajectory(trajectory)
    assert np.isclose(online_trj.calculate(pipeline, data), lindemannindex)


@pytest.mark.parametrize(
    ("trajectory", "lindemannindex"),
    [
        (
            "tests/test_example/459_01.lammpstrj",
            0.025923892565654555,
        ),
        (
            "tests/test_example/459_02.lammpstrj",
            0.026426709832984754,
        ),
    ],
)
def test_frames(trajectory, lindemannindex):
    """Example test with parametrization."""
    frame = read.frames(trajectory)
    test_array = per_frames.calculate(frame)
    assert np.isclose(test_array[-1], lindemannindex)


@pytest.mark.parametrize(
    ("trajectory"),
    [
        ("tests/test_example/459_01.lammpstrj"),
        ("tests/test_example/459_02.lammpstrj"),
    ],
)
def test_frames_middle(trajectory):
    """Example test with parametrization."""
    frame = read.frames(trajectory)
    test_array = per_frames.calculate(frame)
    linde_from_200_trj = per_trj.calculate(frame[0:201])
    linde_at_200_per_frames = test_array[200]
    assert np.isclose(linde_at_200_per_frames, linde_from_200_trj)


@pytest.mark.parametrize(
    ("trajectory", "lindemannindex"),
    [
        (
            "tests/test_example/459_01.lammpstrj",
            0.025923892565654555,
        ),
        (
            "tests/test_example/459_02.lammpstrj",
            0.026426709832984754,
        ),
    ],
)
def test_online_frames(trajectory, lindemannindex):
    """Example test with parametrization."""
    pipeline, data = read.trajectory(trajectory)
    test_array = online_frames.calculate(pipeline, data)
    assert np.isclose(test_array[-1], lindemannindex)


@pytest.mark.parametrize(
    ("trajectory"),
    [
        ("tests/test_example/459_01.lammpstrj"),
        ("tests/test_example/459_02.lammpstrj"),
    ],
)
def test_online_frames_middle(trajectory):
    """Example test with parametrization."""
    pipeline, data = read.trajectory(trajectory)
    test_array = online_frames.calculate(pipeline, data)
    linde_from_200_trj = online_trj.calculate(pipeline, data, 201)
    linde_at_200_per_frames = test_array[200]
    assert np.isclose(linde_at_200_per_frames, linde_from_200_trj)


@pytest.mark.parametrize(
    ("trajectory", "lindemannindex"),
    [
        (
            "tests/test_example/459_01.lammpstrj",
            0.025923892565654555,
        ),
        (
            "tests/test_example/459_02.lammpstrj",
            0.026426709832984754,
        ),
    ],
)
def test_atoms(trajectory, lindemannindex):
    """Example test with parametrization."""
    frame = read.frames(trajectory)
    assert np.isclose(np.mean(per_atoms.calculate(frame)[-1]), lindemannindex)


@pytest.mark.parametrize(
    ("trajectory"),
    [
        ("tests/test_example/459_01.lammpstrj"),
        ("tests/test_example/459_02.lammpstrj"),
    ],
)
def test_atoms_middle(trajectory):
    """Example test with parametrization."""
    frame = read.frames(trajectory)
    per_atoms_frame_at_200 = per_atoms.calculate(frame)[200]
    lindeman_at_200_frame = np.mean(per_atoms_frame_at_200)
    lindeman_from_200_trj = per_trj.calculate(frame[0:201])
    assert np.isclose(lindeman_at_200_frame, lindeman_from_200_trj)
