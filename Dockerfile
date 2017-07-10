FROM ubuntu

RUN apt-get update 

# note what we need to compile neci on a clean ubuntu docker image 

# and i think i should use python2.7-slim, as we need python for the templating
# cmake 
# gfortran
# mpi somehow
# lapack + blas 
# hdf5 maybe.. 
# MPI definetly! especially for the test_suite i need MPI, since it tests it MPI

RUN apt-get install -y cmake 
RUN apt-get install -y gfortran 
RUN apt-get install -y libblas-dev
RUN apt-get install -y liblapack-dev
RUN apt-get install -y build-essential
RUN apt-get install -y python-minimal
RUN apt-get install -y openmpi-bin
RUN apt-get install -y libopenmpi-dev
RUN apt-get install -y git

RUN git clone https://github.com/jsspencer/testcode /testcode

RUN export PYTHONPATH=/testcode/lib:${PYTHONPATH}

# i also have to setup testcode.py and export the correct environment variables
# can i build that from scratch or can i make a docker image or smth.. 
