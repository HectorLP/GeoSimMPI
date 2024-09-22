# GeoSimMPI
This is a modified version of GEOSIM, which generates multi-scale fluvial hetereogeneious facies. This code does not need the same number of cores as the number of compound bars generated. To compile these codes, the Intel fortran and openmpi compilers are both needed. These files need the Intel compiler:

1. cbloc.f90
2. combinefiles.f90
3. IndToVist.f90

These files need the openmpi compiler:

1. main.f90
2. merge.f90

The oder of running the executable files is (with the environment they need):

1. cbloc (Intel compiler environment)
2. main (GCC-OpenMPI environment)
3. merge (GCC-OpenMPI environment)
4. combinefiles (Intel compiler environment)
5. IndToVist (Intel compiler environment)

When the executable files of main and merge, the command is:
mpirun -np N main
mpirun -np N merge
N - replaces N with the number of cores you are going to use.

# The functions of each executable file
1. cbloc: Generate the compound bars in the designated domain and their coordinates and geometric parameters. The locations of cross bar channel fills are also created.
2. main: Generate the unit bars to fill each compound bar and their coordinates and geometric parameters. The polyhydra for cross stratasets are also created in this process.
3. merge: The full channel belt model is stored with all details. The geometrical model files for the sampled area are generated and the sampling resolution is defined by users.
4. combinefiles: All the geometrical mode files are integrated to generate the file including the strata information.
5. IndToVisit: Generate the data file for the visualization in VISIT or PARAVIEW

# Input files and parameters for each executable file
The files listed below do not include those files having the file name in the format: *.out.*
1. cbloc: CBLOC.dat
2. main: UBLOC.dat, UBPLANE.dat, TSLOC.dat
3. merge: indunits.out, sample.dat
4. combinefiles: number of compoun bars
5. IndToVist: indicator.out
To visualize the final result, one more file is needed: merge.bov. User needs to input the file name for the visualization, resolution of the model, and, number of grids.



