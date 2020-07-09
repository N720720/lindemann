# lindemann

<div align="center">

[![Build status](https://github.com/N720720/lindemann/workflows/build/badge.svg?branch=master&event=push)](https://github.com/N720720/lindemann/actions?query=workflow%3Abuild)
[![Python Version](https://img.shields.io/pypi/pyversions/lindemann.svg)](https://pypi.org/project/lindemann/)
[![Dependencies Status](https://img.shields.io/badge/dependencies-up%20to%20date-brightgreen.svg)](https://github.com/N720720/lindemann/pulls?utf8=%E2%9C%93&q=is%3Apr%20author%3Aapp%2Fdependabot)

[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Security: bandit](https://img.shields.io/badge/security-bandit-green.svg)](https://github.com/PyCQA/bandit)  
[![Pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white)](https://github.com/N720720/lindemann/blob/master/.pre-commit-config.yaml)
[![Semantic Versions](https://img.shields.io/badge/%F0%9F%9A%80-semantic%20versions-informational.svg)](https://github.com/N720720/lindemann/releases)
[![License](https://img.shields.io/github/license/N720720/lindemann)](https://github.com/N720720/lindemann/blob/master/LICENSE)

![](459_Atoms_brass.gif)

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

## Usage

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

## Demo

Basic usage to calculate the Lindemann Index:

![](linde_t.gif)

The package has a plotting feature. It will show the a plot Lindemann index vs. the frames. If the trajectory file is a temperature ramp, it is possible to determine the phasetransition.

![](linde_p_new.gif)

Usage of the of the lammpstrj file output feature to save the progression for each atom per frame into a lammps trajectory file. Afterwards the trajectory can be viewed with ovito for example, here the lindemann progression was used for the ovito color coding feature.

![](demo_lammps_ovito.gif)

## Background

A key problem with the measurement of the melting point of nanoparticles is that with decreasing size of a given nanoparticle the phase transition, defined as the temperature of a sudden change in the enthalpy, becomes less pronounced. This is caused by the surface effect: for a given cluster the surface area is larger compared to a bulk structure of the same size. Melting does not take place all at once, but is a longer melt transition and no longer really a melting point. 

The Lindemann index, stated in the following equation presents a solution for this problem. It describes the root-mean-square (rms) fluctuation of the bonds or interatomic distance in the system over time (or temperature, if the temperature of the system changes as the simulation progresses). The Lindemann index is a more robust method to determine the melting point of nanoparticles, opposed to the enthalpy. Accordingly the Lindemann index is often considered, when the melting point of nano-particles is of interest. The index is defined as, 

<a href="https://www.codecogs.com/eqnedit.php?latex={\langle&space;q_{i}&space;\rangle_{\text{atoms}}={\frac&space;{1}{N(N-1)}}\sum&space;_{j\neq&space;i}{\frac&space;{\sqrt&space;{\langle&space;r_{ij}^{2}\rangle&space;-\langle&space;r_{ij}\rangle&space;^{2}}}{\langle&space;r_{ij}\rangle&space;}}}&space;~." target="_blank"><img src="https://latex.codecogs.com/gif.latex?{\langle&space;q_{i}&space;\rangle_{\text{atoms}}={\frac&space;{1}{N(N-1)}}\sum&space;_{j\neq&space;i}{\frac&space;{\sqrt&space;{\langle&space;r_{ij}^{2}\rangle&space;-\langle&space;r_{ij}\rangle&space;^{2}}}{\langle&space;r_{ij}\rangle&space;}}}&space;~." title="{\langle q_{i} \rangle_{\text{atoms}}={\frac {1}{N(N-1)}}\sum _{j\neq i}{\frac {\sqrt {\langle r_{ij}^{2}\rangle -\langle r_{ij}\rangle ^{2}}}{\langle r_{ij}\rangle }}} ~." /></a>

where <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\fn_phv&space;\large&space;N" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;\fn_phv&space;\large&space;N" title="\large N" /></a> is the number of atoms in the nano particle. <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\fn_phv&space;\large&space;r_{ij}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;\fn_phv&space;\large&space;r_{ij}" title="\large r_{ij}" /></a> is the distance between atom i and atom j. The brackets <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\fn_phv&space;\large&space;\langle~\rangle" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;\fn_phv&space;\large&space;\langle~\rangle" title="\large \langle~\rangle" /></a> representing a time or temperature average. The rms of the bond fluctuation is considerable lower for a solid than for a liquid, due to restrained degrees of freedom. In a solid, the atoms hold on to their position and only fluctuate around their equilibrium positions. During the melting process the atoms become more mobile and are able to leave their original position. The translation movement of atoms is magnitudes larger than that of the bond fluctuations of a solid. As a result, the Lindemann index rises dramatically at the melting point and therefore gives a suitable observable to determine the transition phase. In effect, the Lindemann index measures a sort of average difussion coefficient for the atoms in the system.

A key problem with much of the literature regarding the Lindemann index, is that there is a uncertainty of where to define the phase transition within a Lindemann plot. On the grounds that the melting point is a macroscopic definition for bulk structures. But here nanoparticles differ: because of their relative small size, compared to bulk structures, melting can occur fist on their relative large surface, compared to the volume they obtain and followed by the melting of core of the particle. Therefore a temperature range, rather then a melting point, is observed, as stated by Neyts in their work.    

## 🛡 License

[![License](https://img.shields.io/github/license/N720720/lindemann)](https://github.com/N720720/lindemann/blob/master/LICENSE)

This project is licensed under the terms of the `GPLv3` license. See [LICENSE](https://github.com/N720720/lindemann/blob/master/LICENSE) for more details.

## 📃 Citation

```
@misc{lindemann,
  author = {Sebastian Thurm},
  title = {lindemann is a python package to calculate the Lindemann index  of a lammps trajectory as well as the progression of the Lindemann index per frame of temperature ramps  for phase transition analysis.},
  year = {2020},
  publisher = {GitHub},
  journal = {GitHub repository},
  howpublished = {\url{https://github.com/N720720/lindemann}}
}
```

## Credits

The Lindemann index is introduced in the following paper,\
F. A. Lindemann, *Zeitschrift für Phys.* **1910**, *11*, 609–612.\
This package is based on the work from [`ybyygu`](https://github.com/ybyygu/lindemann-index)
and [`whashi44`](https://github.com/whashi44/lindemann) on calculating the Lindemann index.\
The visualisations in this Readme are made with [`Ovito`](https://www.ovito.org/).\
A. Stukowski, *Model. Simul. Mater. Sci. Eng.* **2010**, *18*, 15012. [`link`](https://iopscience.iop.org/article/10.1088/0965-0393/18/1/015012).\
This project was generated with [`python-package-template`](https://github.com/TezRomacH/python-package-template).



