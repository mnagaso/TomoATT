FROM ubuntu:latest AS spython-base
ADD ../../TomoATT /opt/TomoATT
LABEL Maintainer MasaruNagaso
LABEL Version v1.0

ENV PKG_DIR=/opt/ompi
ENV SINGULARITY_PKG_DIR=$PKG_DIR
ENV SINGULARITYENV_APPEND_PATH=$PKG_DIR/bin
ENV SINGULAIRTYENV_APPEND_LD_LIBRARY_PATH=$PKG_DIR/lib
ENV PATH=$PKG_DIR/bin:$PATH
ENV PATH=/opt/TomoATT/build/bin:$PATH

RUN echo "Installing required packages..."
RUN apt-get update && apt-get install -y wget git bash gcc gfortran g++ make file cmake
RUN export PKG_DIR=/opt/packages
RUN export NPROC=4

RUN echo "Installing Open MPI"
RUN export OMPI_VERSOIN_MAJOR=3.0
RUN export OMPI_VERSION_FULL=3.0.1
RUN export OMPI_URL="https://download.open-mpi.org/release/open-mpi/v$OMPI_VERSOIN_MAJOR/openmpi-$OMPI_VERSION_FULL.tar.gz"
RUN mkdir -p /tmp/ompi
RUN mkdir -p /opt
RUN cd /tmp/ompi && wget -O openmpi-$OMPI_VERSION_FULL.tar.gz $OMPI_URL && tar -xvf openmpi-$OMPI_VERSION_FULL.tar.gz
RUN cd /tmp/ompi/openmpi-$OMPI_VERSION_FULL && ./configure --prefix=$PKG_DIR && make -j$NPROC && make install
RUN export PATH=$PKG_DIR/bin:$PATH
RUN export LD_LIBRARY_PATH=$PKG_DIR/lib:$LD_LIBRARY_PATH
RUN export MANPATH=$PKG_DIR/share/man:$MANPATH
RUN echo "openmpi installed."

RUN echo "Compiling hdf5"
RUN export HDF5_VERSION_MAJOR=1.13
RUN export HDF5_VERSION_FULL=1.13.2
RUN export HDF5_URL="https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-$HDF5_VERSION_MAJOR/hdf5-$HDF5_VERSION_FULL/src/hdf5-$HDF5_VERSION_FULL.tar.gz"
RUN mkdir -p /tmp/hdf5
RUN cd /tmp/hdf5 && wget -O hdf5-$HDF5_VERSION_FULL.tar.gz $HDF5_URL && tar -xvf hdf5-$HDF5_VERSION_FULL.tar.gz && cd hdf5-$HDF5_VERSION_FULL
RUN tar -xvf hdf5-$HDF5_VERSION_FULL.tar.gz && cd hdf5-$HDF5_VERSION_FULL
RUN CC=mpicc CXX=mpic++ ./configure --enable-parallel --enable-unsupported --enable-shared --enable-cxx --prefix=$PKG_DIR && make -j$NPROC && make install
RUN echo "hdf5 installed."

RUN cd /opt/TomoATT
RUN mkdir -p build && cd build && rm -rf ./*
RUN cmake .. -DCMAKE_PREFIX_PATH=$PKG_DIR
RUN make -j$NPROC
RUN echo "TomoATT installed."