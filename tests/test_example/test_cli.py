from typer.testing import CliRunner

import lindemann
from lindemann.main import app  # type: ignore

runner = CliRunner()


def test_p_flag():
    result = runner.invoke(app, ["tests/test_example/459_02.lammpstrj", "-p"])
    assert result.exit_code == 0
    assert "Saved file as:" in result.stdout


def single_process_and_multiprocess(trajectory, flag, result_str):
    trajectory.append(flag)
    result = runner.invoke(app, trajectory)
    assert result.exit_code == 0
    assert result_str in result.stdout


def test_all_flags_multiprocess():
    trajectory = ["tests/test_example/459_02.lammpstrj", "tests/test_example/459_01.lammpstrj"]
    result_str = "multiprocessing is implemented only for the -t flag"
    for flag in ["-f", "-of", "-a", "-oa", "-p", "-ti", "-m"]:
        single_process_and_multiprocess(trajectory, flag, result_str)


def test_a_p_flags():
    trajectory = ["tests/test_example/459_02.lammpstrj"]
    result_str = "lindemann_index_per"
    for flag in [
        "-a",
        "-oa",
        "-f",
        "-of",
        "-p",
    ]:
        single_process_and_multiprocess(trajectory, flag, result_str)


def test_t_flag():
    flag = "-t"
    res_str = "lindemann index for the Trajectory: 0.026426"
    trajectory = ["tests/test_example/459_02.lammpstrj"]
    single_process_and_multiprocess(trajectory, flag, res_str)


def test_m_flag():
    flag = "-m"
    res_str = "Memory use: \nFlag -t (per_trj) will use 0.0034 GB\nFlag -ot (per_trj) will use 0.0008 GB\nFlag -f (per_frames) will use 0.0034 GB\nFlag -of (per_frames) will use 0.0008 GB\nFlag -a (per_atoms) will use 0.005 GB\nFlag -oa (per_atoms) will use 0.0024 GB\n"
    trajectory = ["tests/test_example/459_02.lammpstrj"]
    single_process_and_multiprocess(trajectory, flag, res_str)


def test_none_flag():
    result = runner.invoke(app, ["tests/test_example/459_02.lammpstrj"])
    assert result.exit_code == 0
    assert "lindemann index for the Trajectory: 0.026426" in result.stdout


def test_none_flag_multi():
    result = runner.invoke(
        app, ["tests/test_example/459_01.lammpstrj", "tests/test_example/459_02.lammpstrj"]
    )
    assert result.exit_code == 0
    assert "0.025923" in result.stdout
    assert "0.026426" in result.stdout
