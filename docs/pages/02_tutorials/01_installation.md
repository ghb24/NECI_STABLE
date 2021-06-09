title: NECI Installation Guide
author: Philip Haupt
---

@note
This is modified from instructions found on the internal wiki, written by Kai and Oskar.
@endnote

# NECI Installation

Installation of NECI requires[^allogin]

-   git,
-   cmake 3.12 or newer,
-   BLAS/LAPACK,
-   MPI 3,
-   a Fortran compiler supporting Fortran 2003 features, and
-   optionally HDF5 (recommended).

There are two git repositories for NECI. The stable release is available publicly on github, [here](https://github.com/ghb24/NECI_STABLE).
The developer version is on a private repository on bitbucket, for which you need to be invited to see.
The Github version is not updated as frequently, so if you wish to use the latest methods, contact the developers.

To get started, we must first clone the repository, with
```bash
https://github.com/ghb24/NECI_STABLE.git
```
for the public release, or
```bash
git clone https://<username>@bitbucket.org/neci_developers/neci.git
```
for the developer release (replace `<username>` with your bitbucket username).

Next, make a subdirectory (let's call it "build"), then run cmake and then make:
```bash
mkdir build
cd build
cmake -DENABLE_HDF5=ON -DENABLE_BUILD_HDF5=ON ../
make -j
```
@note
if you are making without HDF5, then remove the options `-DENABLE_HDF5=ON -DENABLE_BUILD_HDF5=ON` from the third line.
@endnote

## HDF5

You may also need to build HDF5 yourself as a shared library, for example to ensure the same compilers for HDF5 and your programs.

In this case, download and extract the program from the HDF5 group website, then build with (you may replace the compilers if you wish, for example use `FC=mpifort` instead of `FC=mpiifort`)
```bash
cd your_build_directory
export HDF5_SRC= # The HDF5 source
export HDF5_ROOT= # Where it should be installed
FC=mpiifort F9X=mpiifort CC=mpicc $HDF5_SRC/configure --prefix=$HDF5_ROOT --enable-fortran --enable-fortran2003 --enable-parallel
make
make install
```
where you define `HDF5_SRC` and `HDF5_ROOT` appropriately. Then, before running `cmake` for NECI, run

```bash
export HDF5_ROOT= # where HDF5 was installed in the previous step
```

and proceed with the NECI installation as before, 
```bash
mkdir build
cd build
cmake -DENABLE_HDF5=ON ../
make -j
```


[^allogin]: If you are on a cluster, you may need to run a command similar to `module load ifort mpi.intel`.

///Footnotes Go Here///
