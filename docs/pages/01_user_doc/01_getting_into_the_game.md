---
title: Installation
---


## Getting the Code


There are two git repositories for NECI. The stable release is available publicly on github, [here](https://github.com/ghb24/NECI_STABLE).
The developer version is on a private repository on bitbucket, for which you need to be invited to see.
The Github version is not updated as frequently, so if you wish to use the latest methods, contact the developers.

### BitBucket Repository Access

To gain access, contact one of the repository administrators [Oskar Weser
([o.weser@fkf.mpg.de](mailto:o.weser@fkf.mpg.de)) and Werner Dobrautz
([w.dobrautz@fkf.mpg.de](mailto:w.dobrautz@fkf.mpg.de))] who will invite you. If you
already have a bitbucket account let the repository administrators know
the email address associated with your account.

You will receive an invitation email. Please accept this invitation, and
create a bitbucket account as prompted if necessary.

To gain access to the NECI repository, an ssh key is required. This can
be generated on any linux machine using the command\footnote{`ssh-keygen` can also generate DSA keys. Some ssh clients and servers will reject DSA keys longer than 1024 bits, and 1024 bits is	currently on the margin of being crackable. As such 2048 bit RSA keys are preferred. Top secret this code is. Probably. Apart from the master branch which hosted for all on github. And in molpro.	And anyone that wants it obviously.}

```bash
ssh-keygen -t rsa -b 2048
```

This will create a private (`~/.ssh/id_rsa`) and a public key file
(`~/.ssh/id_rsa.pub`).

The private key must be kept private. On the bitbucket homepage, go to
account settings (accessible from the top-right of the main page), and
navigate to “SSH keys”. Click “Add key” and add the contents of the
public key. This will give you access to the repository.

You can now clone the code into a new directory using the command

```bash
git clone git@bitbucket.org:neci_developers/neci.git [target_dir]
```

## NECI Installation

@note
This is modified from instructions found on the internal wiki, written by Kai and Oskar.
@endnote

Installation of NECI requires\footnote{If you are on a cluster, you may need to run a command similar to `module load ifort mpi.intel`.}

-   git,
-   cmake 3.12 or newer,
-   BLAS/LAPACK,
-   MPI 3,
-   a Fortran compiler supporting Fortran 2003 features, and
-   optionally HDF5 (recommended).

To get started, we must first clone the repository, with
```bash
https://github.com/ghb24/NECI_STABLE.git
```
for the public release, or
```bash
git clone https://<username>@bitbucket.org/neci_developers/neci.git
```
for the developer release (replace `<username>` with your bitbucket username).
The directory of this repository will be referred to as "`your_code_directory`".

Next, create a directory (let's call it "`your_build_directory`"), then run cmake and then make:
```bash
mkdir build
cd build
cmake -DENABLE_HDF5=ON -DENABLE_BUILD_HDF5=ON ${your_code_directory}
make -j
```
@note
If you are making without HDF5, then remove the options `-DENABLE_HDF5=ON -DENABLE_BUILD_HDF5=ON` from the third line.
@endnote

## HDF5

You may also need to build HDF5 yourself as a shared library, which speeds up the compilation process, since HDF5 does not have to be rebuilt for every new project.
@note
HDF5 should be built with the same set of compilers as NECI.
@endnote

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
cmake -DENABLE_HDF5=ON ${your_code_directory}
make -j
```

## Documentation

If you wish to generate documentation for NECI, based on [Ford](https://github.com/Fortran-FOSS-Programmers/ford),
then when you build with CMake, you must also include the flag `-DENABLE_DOC=ON`, for example
```
mkdir build
cd build
cmake -DENABLE_DOC=ON ${your_code_directory}
make -j doc
```

Requirements to produce the docs are:

- pandoc,
- latexmk,
- pdflatex,
- biber, and
- a working internet connection (only for the first time, in order to get a custom-build of Ford).

This will not only generate this documentation in the form of a PDF, but also as a website,
having in addition to the information in the PDF also automatically generated documentation from comments in the source files.
