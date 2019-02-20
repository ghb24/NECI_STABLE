FROM opensuse:42.3

RUN zypper update 

# note what we need to compile neci on a clean ubuntu docker image 

# and i think i should use python2.7-slim, as we need python for the templating
# cmake 
# gfortran
# mpi somehow
# lapack + blas 
# hdf5 maybe.. 
# MPI definetly! especially for the test_suite i need MPI, since it tests it MPI

RUN zypper install -y cmake 
RUN zypper install -y gcc-fortran 
RUN zypper install -y blas-devel
RUN zypper install -y lapack-devel
RUN zypper install -y gcc-c++
# i still have to figure out what i need instead of this:
#RUN apt-get install -y build-essential
RUN zypper install -y python
RUN zypper install -y openmpi
RUN zypper install -y openmpi-devel
# but somehow zypper does not add this to the path.. so do it manually
ENV PATH=/usr/lib64/mpi/gcc/openmpi/bin:${PATH}
ENV INCLUDE=/usr/lib64/mpi/gcc/openmpi/include:/usr/lib64/mpi/gcc/openmpi/lib64:${INCLUDE}
ENV LD_LIBRARY_PATH=/usr/lib64/mpi/gcc/openmpi/lib64:${LD_LIBRARY_PATH}

RUN zypper install -y git
# already included in cmake package:
#RUN apt-get install -y cmake-curses-gui
# to enable hdf5 compilation with cmake we also need: 
RUN zypper install -y zlib zlib-devel

RUN git clone https://github.com/jsspencer/testcode /testcode

RUN export PYTHONPATH=/testcode/lib:${PYTHONPATH}

# i also have to setup testcode.py and export the correct environment variables
# can i build that from scratch or can i make a docker image or smth.. 
