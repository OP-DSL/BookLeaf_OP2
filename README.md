# BookLeaf


## Introduction

Bookleaf is an unstructured Lagrangian Hydro mini-app.

Four input decks are provided: Sod, Sedov, Saltzmann and Noh.

Current BookLeaf_ref as V1.1, plus:
* parallel ALE capability
* Modified makefile.intel
* Input endtimes modified to match test runs


Wiki: https://github.com/UK-MAC/BookLeaf/wiki


## BookLeaf Build Procedure

BookLeaf uses OP2-Common's current Fortran runtime and translator-v2. The
root `Makefile` is the only build file. Translator output is written to
`generated/bookleaf/`; it is a build product.

Initialise the configured OP2 environment before building:

```bash
source ../OP2-Common/scripts/source_gnuz
make
```

Individual targets are:

- `bookleaf_seq`: direct developer build, without translation.
- `bookleaf_genseq`: translator-v2 sequential backend.
- `bookleaf_openmp`: translator-v2 OpenMP backend.
- `bookleaf_c_cuda`, `bookleaf_c_hip`: translator-v2 C/CUDA or C/HIP backend.
- `bookleaf_mpi_seq`, `bookleaf_mpi_genseq`, `bookleaf_mpi_openmp`,
  `bookleaf_mpi_c_cuda`, `bookleaf_mpi_c_hip`: corresponding MPI builds.

`make` includes only configured variants. Every named target can also be
built explicitly; accelerator targets require the relevant compiler and OP2
library to have been configured, pointed by OP2_INSTALL_PATH.

## Generating a mesh

The mesh generator is built directly from BookLeaf's original mesh-construction
sources.

```bash
make bookleaf_meshgen
./bookleaf_meshgen FILE=input/noh
```

This writes `mesh_hdf5` in the current directory. Generate it in a separate
working directory if an existing mesh file must be retained.

## Running the Code

The application reads `mesh_hdf5` from its current working directory.

```bash
./bookleaf_seq FILE=input/noh
./bookleaf_genseq FILE=input/noh
OMP_NUM_THREADS=4 ./bookleaf_openmp FILE=input/noh
./bookleaf_c_cuda FILE=input/noh
./bookleaf_c_hip FILE=input/noh

mpirun -np 4 ./bookleaf_mpi_seq FILE=input/noh
mpirun -np 4 ./bookleaf_mpi_genseq FILE=input/noh
OMP_NUM_THREADS=8 mpirun -np 4 --map-by ppr:4:node:PE=8 --bind-to core ./bookleaf_mpi_openmp FILE=input/noh
mpirun -np 4 ./bookleaf_mpi_c_cuda FILE=input/noh
mpirun -np 4 ./bookleaf_mpi_c_hip FILE=input/noh
```

## Version History

BookLeaf_ref - As V1.2

V1.2   - Adds in parallel ALE. Plus:
* Modified makefile.intel for Xeon vectorisation at OPT level and PHI=1 option to build for Xeon Phi, -qopt-report=3 no longer default flag
* Makefile help has PHI option added plus version updated to v1.1
* End times for sod and sedov test cases changed to reflect test cases run

V1.1   - Adds in mesh partitioning. Parallel running now available.

V1.0   - Initial version. Contains MPI comms, but only serial meshes can be contructed.
