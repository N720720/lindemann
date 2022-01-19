# `lindemann`

lindemann is a python package to calculate the Lindemann index  of a lammps trajectory as well as the progression of the Lindemann index per frame of temperature ramps  for phase transition analysis.

**Usage**:

```console
$ lindemann [OPTIONS] TRJFILE...
```

**Arguments**:

* `TRJFILE...`: The trajectory file(s). If no other option is selected, the lindemann index is calculated for the trajectory. Equivalent to the -t option. If you pass more than one trajectory they will be calculated in parallel. Only works with no flag or -t flag.   [required]

**Options**:

* `-t`: Calculates the Lindemann-Index for the Trajectory file(s)  [default: False]
* `-f`: Calculates the Lindemann-Index for each frame.  [default: False]
* `-a`: Calculates the Lindemann-Index for each atom for each frame.  [default: False]
* `-p`: Returns a plot Lindemann-Index vs. Frame.  [default: False]
* `-l`: Saves the individual Lindemann-Index of each Atom in a lammpstrj, so it can be viewed in Ovito.  [default: False]
* `-v, --version`: Prints the version of the lindemann package.
* `-ti, -timeit`: Uses timeit module to show running time  [default: False]
* `-m, -mem_use`: Calculates the memory use. Run it before you use any of the cli functionality despite the -t flag  [default: False]
* `--help`: Show this message and exit.
