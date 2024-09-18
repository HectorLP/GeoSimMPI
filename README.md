# GeoSimMPI
This is a modified version of GEOSIM, which generates multi-scale fluvial hetereogeneious facies. This code does not need the same number of cores as the number of compound bars generated. To compile these codes, the Intel fortran and openmpi compilers are both needed. These files need the Intel compiler:

1. cbloc.f90
2. combinefiles.f90
3. IndToVist.f90

These files need the openmpi compiler:

1. main.f90
2. merge.f90

The oder of running the executable files is:

1. cbloc
2. main
3. merge
4. combinefiles
5. IndToVist (if you need to visualize the structure in the Visit)
