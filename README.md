# lindemann

<div align="center">

[![Build status](https://github.com/N720720/lindemann/workflows/build/badge.svg?branch=master&event=push)](https://github.com/N720720/lindemann/actions?query=workflow%3Abuild)
[![Python Version](https://img.shields.io/pypi/pyversions/lindemann.svg)](https://pypi.org/project/lindemann/)
[![Dependencies Status](https://img.shields.io/badge/dependencies-up%20to%20date-brightgreen.svg)](https://github.com/N720720/lindemann/pulls?utf8=%E2%9C%93&q=is%3Apr%20author%3Aapp%2Fdependabot)

[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Security: bandit](https://img.shields.io/badge/security-bandit-green.svg)](https://github.com/PyCQA/bandit)
[![Pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white)](https://github.com/N720720/lindemann/blob/master/.pre-commit-config.yaml)
[![Semantic Versions](https://img.shields.io/badge/%F0%9F%9A%80-semantic%20versions-informational.svg)](https://github.com/N720720/lindemann/releases)
[![License](https://img.shields.io/badge/License-GPL--3.0 License-yellow.svg)](https://github.com/N720720/lindemann/blob/master/LICENSE)

lindemann is a python package to calculate the Lindemann index  of a lammps trajectory as well as the progression of the Lindemann index per frame of temperature ramps  for phase transition analysis.
</div>


## Installation

```bash
pip install lindemann
```

or install with `Poetry`

```bash
poetry add lindemann
```

**Usage**:

```console
$ lindemann [OPTIONS] TRJFILE
```
or if installed with `Poetry`:

```bash
poetry run [OPTIONS] TRJFILE
```
**Options**:

* `-t`: Calculates the Lindemann-Index for the Trajectory file
* `-f`: Calculates the Lindemann-Index for each frame.
* `-a`: Calculates the Lindemann-Index for each atom for each frame.
* `-p`: Returns a plot Lindemann-Index vs. Frame.
* `-l`: Saves the individual Lindemann-Index of each Atom in a lammpstrj, so it can be viewed in Ovito.
* `-v, --version`: Prints the version of the lindemann package.
* `-ti, -timeit`: Uses timeit module to show running time
* `--help`: Show this message and exit.

## 🛡 License

[![License](https://img.shields.io/github/license/N720720/lindemann)](https://github.com/N720720/lindemann/blob/master/LICENSE)

This project is licensed under the terms of the `MIT` license. See [LICENSE](https://github.com/N720720/lindemann/blob/master/LICENSE) for more details.

## 📃 Citation

```
@misc{lindemann,
  author = {N720720},
  title = {lindemann is a python package to calculate the Lindemann index  of a lammps trajectory as well as the progression of the Lindemann index per frame of temperature ramps  for phase transition analysis.},
  year = {2020},
  publisher = {GitHub},
  journal = {GitHub repository},
  howpublished = {\url{https://github.com/N720720/lindemann}}
}
```

## Credits

This project was generated with [`python-package-template`](https://github.com/TezRomacH/python-package-template).
