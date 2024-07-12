#!/usr/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=04:00:00
#SBATCH --mail-type=none
#SBATCH -p sapphire
#SBATCH -A ukaea-ap001-CPU
#SBATCH --cpus-per-task=1
#SBATCH -o out_%j_%A_%a
#SBATCH --exclusive

. /etc/profile.d/module.sh
module purge
module load rhel8/slurm dot rhel8/global
module load gcc/11
module load python/3.11.0-icl
#Did build with this but seemed a bit flaky
#but could try this intead of building an MPI
# module load openmpi/gcc/9.3/4.0.4

#run sbatch <script-name> from a directory where you want firedrake
#after it runs it should have
#a usr directory where openmpi is installed
#a firedrake directory
#a env.sh file that can be sourced to get the build environment back

ROOTDIR=${SLURM_SUBMIT_DIR}

#Build openMPI + download if it isn't here already
if [ ! -d openmpi-4.1.6 ]; then
  if [ ! -f openmpi-4.1.6.tar.gz ]; then
    wget https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-4.1.6.tar.gz \
        || { echo " Failed to download openmpi" ; exit 1 ; }
    tar -xf openmpi-4.1.6.tar.gz
  fi
fi

cd openmpi-4.1.6
CC=gcc CXX=g++ FC=gfortran ./configure --with-pmix=internal --prefix=${ROOTDIR}/usr \
        || { echo " Failed to configure openmpi" ; exit 1 ; }

make -j 20  || { echo " Failed to build openmpi" ; exit 1 ; }
make install || { echo " Failed to install openmpi" ; exit 1 ; }


export PATH=${ROOTDIR}/usr/bin:${PATH}
export LD_LIBRARY_PATH=${ROOTDIR}/usr/lib:${LD_LIBRARY_PATH}
export LIBRARY_PATH=${ROOTDIR}/usr/lib:${LIBRARY_PATH}
export CPATH=${ROOTDIR}/usr/include:${CPATH}
export CPLUS_INCLUDE_PATH=${ROOTDIR}/usr/include:${CPLUS_INCLUDE_PATH}
export C_INCLUDE_PATH=${ROOTDIR}/usr/include:${C_INCLUDE_PATH}


cd ${ROOTDIR}

export CC=mpicc
export CXX=mpicxx
export FC=mpif90
export F77=mpif70

#Have in the past had trouble using -march=native for some mystery reason
#it doesn't like avx512 though it might work?
#Also not 100% sure these are passing through to petsc
export COPTFLAGS='-O3 -mavx2' 
export CXXOPTFLAGS='-O3 -mavx2' 
export FOPTFLAGS='-O3 -mavx2'

#Get the firedrake script
if [ ! -f firedake-install ]; then 
	curl -O https://raw.githubusercontent.com/firedrakeproject/firedrake/master/scripts/firedrake-install
fi
python3 firedrake-install --venv-name firedrake --disable-ssh --mpiexec=mpirun --mpicc=$CC --mpicxx=$CXX --mpif90=$FC --no-package-manager --install irksome


#Create a file to source to get the same env back
echo "#!/usr/bin/bash" > env.sh
echo ". /etc/profile.d/module.sh" >> env.sh
echo "module purge" >> env.sh
echo "module load rhel8/slurm dot rhel8/global" >> env.sh
echo "module load gcc/11" >> env.sh
echo "module load python/3.11.0-icl" >> env.sh

echo "export PATH=\${ROOTDIR}/usr/bin:\${PATH}" >> env.sh
echo "export LD_LIBRARY_PATH=\${ROOTDIR}/usr/lib:\${LD_LIBRARY_PATH}" >> env.sh
echo "export LIBRARY_PATH=\${ROOTDIR}/usr/lib:\${LIBRARY_PATH}" >> env.sh
echo "export CPATH=\${ROOTDIR}/usr/include:\${CPATH}" >> env.sh
echo "export CPLUS_INCLUDE_PATH=\${ROOTDIR}/usr/include:\${CPLUS_INCLUDE_PATH}" >> env.sh
echo "export C_INCLUDE_PATH=\${ROOTDIR}/usr/include:\${C_INCLUDE_PATH}" >> env.sh