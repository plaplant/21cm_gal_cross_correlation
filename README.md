# 21cm--galaxy Cross-Correlation
This repo contains methods for calculating projected sensitivity for 21cm x
galaxy observations, and a bit of sample data which can be used to derive some
of the results of [La Plante et
al. 2022](https://ui.adsabs.harvard.edu/abs/2022arXiv220509770L/abstract). Most
of the calculation is done in Fortran, though there are portions in Python, and
an accompanying [Jupyter notebook](sensitivity_calc.ipynb) that shows some of
the results. Note that for space reasons, we only include the raw data for our
fiducial set of observing properties. Additional simulation results are
available upon request.

## Compiling the Fortran code
To compile the Fortran code, you must have either `gfortran` or `ifort`
available on your system. Also, you need access to a recent-ish version of HDF5
to read the input data, and write files. No other external dependencies are
required.

Before compiling, you may have to adjust the Makefile to point to your local
HDF5 installation. You should only have to edit the line defining `H5HOME`, as
everything else (compiler, include files, and libraries) are defined relative to
that.

To compile the program, run
```
make
```
This will make the executable `compute_snr.x`. To change any of the observing
parameters or paths to files, edit the values in `compute_snr.f90`.

## Data files
This repo includes some data files that are needed to reproduce the calculation.
They are located in the `data` directory. More information about the files is
available there.
