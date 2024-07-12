# DockerFile for a plain Firedrake suitable for testing Firedrake components and applications

FROM petsc-env

# What is this container?
LABEL name="firedrake", \
    description="A minimal Firedrake installation"

# Set this to 'master' to use current firedrake.
# Set this to a different branch to build from
# there, which can include modifications to the
# install script.
ARG FIREDRAKE_BRANCH=master

LABEL firedrake.branch=${FIREDRAKE_BRANCH}

USER firedrake
WORKDIR /home/firedrake

# Copy the firedrake-install from the correct branch
RUN curl -O https://raw.githubusercontent.com/firedrakeproject/firedrake/${FIREDRAKE_BRANCH}/scripts/firedrake-install
RUN bash -c 'python3 firedrake-install \
    --no-package-manager \
    --disable-ssh \
    --honour-petsc-dir \
    --install irksome \
    --mpicc=$MPICH_DIR/mpicc \
    --mpicxx=$MPICH_DIR/mpicxx \
    --mpif90=$MPICH_DIR/mpif90 \ 
    --mpiexec=$MPICH_DIR/mpiexec \
    --package-branch firedrake ${FIREDRAKE_BRANCH}'

RUN bash -c '. /home/firedrake/firedrake/bin/activate; \
    cd /home/firedrake/firedrake/src/firedrake/; \
    python -m pytest -v tests/regression/ \
    -k "poisson_strong or stokes_mini or dg_advection"; sleep 5'

RUN bash -c '. /home/firedrake/firedrake/bin/activate; \
    firedrake-clean'