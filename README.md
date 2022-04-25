# TopOpt_in_PETSc_Transient
===============
A 3D large-scale TRANSIENT topology optimization code using PETSc
===============

The code (or framework) presented on this page is a fully parallel framework for conducting very large scale transient topology optimziation on structured grids. For more details see www.topopt.dtu.dk/PETSc.

April, 2022, Hansotto Kristiansen & Niels Aage

To clone repository:

        git clone https://github.com/topopt/TopOpt_in_PETSc_Transient.git

NOTE: The code requires PETSc version 3.14.0 or newer ! Also note that the code is not tested against the development branch on git.

This code has been tested on:

    Linux systems including: Ubuntu 18.04, Red hat enterprise linux 8

This code requires the following external software to work:

    PETSc version 3.11.4 or earlier (though never than 3.8.x)
    Requires LAPACK/BLAS
    Requires MPI

Compile following rules in makefile_ref

Normal compilation time of framework, e.g. 4s: "make topopt -j"

Run the base example by typing e.g.: "mpirun -np 4 ./topopt"

Postprocess results using Python 2.6: "bin2vtu #" where # refers to the iteration number

Visualize using ParaView (version 5.7 or earlier)

Paper reference will follow here once final acceptance is obtained
