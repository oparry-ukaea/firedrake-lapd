# DockerFile for an environment into which firedrake can be installed.

# setup the base OS
FROM ubuntu:22.04
USER root
SHELL ["/bin/bash", "-o", "pipefail", "-c"]

# What is this container?
LABEL name="petsc-environment", \
    description="A PETSc build with a minimal set of packages for Firedrake"

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# User inputs: build variables

# compiler flags
ARG BUILD_ARCH="native"
ARG BUILD_COPTFLAGS='-O3 -march=$BUILD_ARCH -mtune=$BUILD_ARCH'
ARG BUILD_CXXOPTFLAGS='-O3 -march=$BUILD_ARCH -mtune=$BUILD_ARCH'
ARG BUILD_FOPTFLAGS='-O3 -march=$BUILD_ARCH -mtune=$BUILD_ARCH'

# how many cores can make use
ARG BUILD_MAKE_NP=8

# debug build?
ARG BUILD_DEBUG=0

# what mpi to use
ARG MPICH_DOWNLOAD_VERSION="4.1.2"
ARG MPICH_DOWNLOAD_DEVICE="ch3:sock"

# extra packages for petsc to install e.g. netcdf, suitesparse, ...
ARG EXTRA_PACKAGES=""

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# essential packages installed on every build
ENV CORE_PACKAGES=" \
    hdf5 \
    hwloc \
    mumps \
    ptscotch \
    scalapack"

LABEL build.debug=${BUILD_DEBUG} \
    build.mpi.version=${MPICH_DOWNLOAD_VERSION} \
    build.mpi.device=${MPICH_DOWNLOAD_DEVICE} \
    build.extra-packages=${EXTRA_PACKAGES}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Update and install required packages for Firedrake
RUN apt-get update \
    && apt-get -y dist-upgrade \
    && DEBIAN_FRONTEND=noninteractive apt-get -y install tzdata \
    && apt-get -y install curl vim docker.io \
    openssh-client build-essential autoconf automake \
    cmake gfortran git libopenblas-serial-dev \
    libtool python3-dev python3-pip python3-tk python3-venv \
    python3-requests zlib1g-dev libboost-dev sudo gmsh \
    bison flex \
    liboce-ocaf-dev \
    swig graphviz \
    libcurl4-openssl-dev libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

# Use a more sane locale
ENV LC_ALL C.UTF-8

# Set up user so that we do not run as root
RUN useradd -m -s /bin/bash -G sudo firedrake && \
    echo "firedrake:docker" | chpasswd && \
    echo "firedrake ALL=(ALL) NOPASSWD: ALL" >> /etc/sudoers && \
    ldconfig

USER firedrake
WORKDIR /home/firedrake

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Fetch PETSc
RUN git clone https://github.com/firedrakeproject/petsc.git

# Build MPICH manually because we don't want PETSc to build it twice
RUN bash -c 'cd petsc; \
    export MPICH_DOWNLOAD_URL=https://github.com/pmodels/mpich/releases/download/v${MPICH_DOWNLOAD_VERSION}/mpich-${MPICH_DOWNLOAD_VERSION}.tar.gz; \
    DOWNLOAD_PACKAGES=" "; \
    for pkg in $CORE_PACKAGES $EXTRA_PACKAGES; \
    do DOWNLOAD_PACKAGES+="--download-$pkg "; done; \
    ./configure \
    --COPTFLAGS=$BUILD_COPTFLAGS \
    --CXXOPTFLAGS=$BUILD_CXXOPTFLAGS \
    --FOPTFLAGS=$BUILD_FOPTFLAGS \
    --with-c2html=0 \
    --with-debugging=$BUILD_DEBUG \
    --with-fortran-bindings=0 \
    --with-make-np=$BUILD_MAKE_NP \
    --with-shared-libraries=1 \
    --with-zlib \
    --download-mpich=$MPICH_DOWNLOAD_URL \
    --download-mpich-device=$MPICH_DOWNLOAD_DEVICE \
    $DOWNLOAD_PACKAGES \
    PETSC_ARCH=packages; \
    mv packages/include/petscconf.h packages/include/old_petscconf.nope;'
# Don't run make here, we only want MPICH and HWLOC
# It is also necessary to move `petscconf.h` so packages isn't treated like a working PETSc

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Build default Firedrake PETSc
RUN bash -c 'cd petsc; \
    export PACKAGES=/home/firedrake/petsc/packages; \
    DOWNLOAD_PACKAGES_DIR=" "; \
    for pkg in $CORE_PACKAGES $EXTRA_PACKAGES; \
    do DOWNLOAD_PACKAGES_DIR+="--with-$pkg-dir=$PACKAGES "; done; \
    ./configure \
    --COPTFLAGS=$BUILD_COPTFLAGS \
    --CXXOPTFLAGS=$BUILD_CXXOPTFLAGS \
    --FOPTFLAGS=$BUILD_FOPTFLAGS \
    --with-c2html=0 \
    --with-debugging=$BUILD_DEBUG \
    --with-fortran-bindings=0 \
    --with-make-np=$BUILD_MAKE_NP \
    --with-shared-libraries=1 \
    --with-bison \
    --with-flex \
    --with-zlib \
    --with-mpi-dir=$PACKAGES \
    $DOWNLOAD_PACKAGES_DIR \
    PETSC_ARCH=default; \
    make PETSC_DIR=/home/firedrake/petsc PETSC_ARCH=default all;'

RUN bash -c 'cd petsc; \
    make PETSC_DIR=/home/firedrake/petsc PETSC_ARCH=default check;'

# Clean up unnecessary files
RUN rm -rf /home/firedrake/petsc/**/externalpackages \
    && rm -rf /home/firedrake/petsc/src/docs \
    && rm -f /home/firedrake/petsc/src/**/tutorials/output/* \
    && rm -f /home/firedrake/petsc/src/**/tests/output/*

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Set some useful environment variables
ENV PETSC_ARCH=default \
    PETSC_DIR=/home/firedrake/petsc \
    MPICH_DIR=/home/firedrake/petsc/packages/bin \
    HDF5_DIR=/home/firedrake/petsc/packages \
    HDF5_MPI=ON \
    OMP_NUM_THREADS=1 \
    OPENBLAS_NUM_THREADS=1