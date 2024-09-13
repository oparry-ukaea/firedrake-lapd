#!/bin/bash

# Slurm job options (name, compute nodes, job time)
#SBATCH --job-name=fd-singularity
#SBATCH --time=04:00:00
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=32
#SBATCH --exclusive
#SBATCH --cpus-per-task=1
#SBATCH --partition=standard
#SBATCH --qos=standard
#SBATCH --account=[SET_ACCOUNT_HERE]

# ===================== Options =====================
# Script filename - **must be in $working_dir!**
script="[SET_SCRIPT_NAME_HERE]"
# Directory containing the singularity image
sif_dir=[SET_IMAGE_DIR_HERE]
# Singularity image name
sif_fname="firedrake.sif"
# Working directory for the run where output will be written (bound to working directory in the container)
working_dir="${SLURM_SUBMIT_DIR}"

# Shouldn't need to change anything below...
# ===================== Modules =====================
module purge
module load load-epcc-module
module load PrgEnv-gnu
module swap cray-mpich cray-mpich-abi
module load cray-dsmml
module load cray-libsci
module load xpmem

# Make srun inherit cpu-per-task settings (if running via salloc)
export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK

# ===================== Environment setup for singularity =====================
# Set the LD_LIBRARY_PATH environment variable within the Singularity container
# to ensure that it used the correct MPI libraries.
export SINGULARITYENV_LD_LIBRARY_PATH="/opt/cray/pe/mpich/8.1.23/ofi/gnu/9.1/lib-abi-mpich:/opt/cray/pe/mpich/8.1.23/gtl/lib:/opt/cray/libfabric/1.12.1.2.2.0.0/lib64:/opt/cray/pe/gcc-libs:/opt/cray/pe/gcc-libs:/opt/cray/pe/lib64:/opt/cray/pe/lib64:/opt/cray/xpmem/default/lib64:/usr/lib64/libibverbs:/usr/lib64:/usr/lib64"

# This makes sure HPE Cray Slingshot interconnect libraries are available
# from inside the container.
export SINGULARITY_BIND="/opt/cray,/var/spool,/opt/cray/pe/mpich/8.1.23/ofi/gnu/9.1/lib-abi-mpich:/opt/cray/pe/mpich/8.1.23/gtl/lib,/etc/host.conf,/etc/libibverbs.d/mlx5.driver,/etc/libnl/classid,/etc/resolv.conf,/opt/cray/libfabric/1.12.1.2.2.0.0/lib64/libfabric.so.1,/opt/cray/pe/gcc-libs/libatomic.so.1,/opt/cray/pe/gcc-libs/libgcc_s.so.1,/opt/cray/pe/gcc-libs/libgfortran.so.5,/opt/cray/pe/gcc-libs/libquadmath.so.0,/opt/cray/pe/lib64/libpals.so.0,/opt/cray/pe/lib64/libpmi2.so.0,/opt/cray/pe/lib64/libpmi.so.0,/opt/cray/xpmem/default/lib64/libxpmem.so.0,/run/munge/munge.socket.2,/usr/lib64/libibverbs/libmlx5-rdmav34.so,/usr/lib64/libibverbs.so.1,/usr/lib64/libkeyutils.so.1,/usr/lib64/liblnetconfig.so.4,/usr/lib64/liblustreapi.so,/usr/lib64/libmunge.so.2,/usr/lib64/libnl-3.so.200,/usr/lib64/libnl-genl-3.so.200,/usr/lib64/libnl-route-3.so.200,/usr/lib64/librdmacm.so.1,/usr/lib64/libyaml-0.so.2"

# Set environment variables inside the Singularity container for firedrake et al.
# Don't multithread
export SINGULARITYENV_OMP_NUM_THREADS=1
# Use the mpi compilers from the firedrake container
export SINGULARITYENV_PYOP2_CC=/home/firedrake/firedrake/bin/mpicc
export SINGULARITYENV_PYOP2_CXX=/home/firedrake/firedrake/bin/mpicxx
# Save caches locally so they persist across `singularity run` calls
export SINGULARITYENV_PYOP2_CACHE_DIR=/home/firedrake/work/.cache/pyop2_cache
export SINGULARITYENV_FIREDRAKE_TSFC_KERNEL_CACHE_DIR=/home/firedrake/work/.cache/tsfc

#==============================================================================
# Check variables have been set
for var in "script" "sif_dir" "sif_fname" "working_dir"; do
     if [ -z "${!var}" ]; then
          echo "'${var}' option wasn't set"
          exit 1
     fi
done
#==============================================================================
# Fudge to force initalisation in serial and bypass lock file issues
echo "Initialising container env"
singularity exec --bind ${working_dir}:/home/firedrake/work --home ${working_dir} "${sif_dir}/${sif_fname}" /home/firedrake/firedrake/bin/python -c "import firedrake"

echo
echo "Running ${script} in ${working_dir} with singularity image ${sif_dir}/${sif_fname}"
# Launch the parallel job
srun --hint=nomultithread --distribution=block:block --chdir=${working_dir}\
     singularity exec \
     --bind ${working_dir}:/home/firedrake/work \
     --home ${working_dir} \
     "${sif_dir}/${sif_fname}" \
     /home/firedrake/firedrake/bin/python \
     "${script}"
