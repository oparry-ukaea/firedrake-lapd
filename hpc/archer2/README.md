## Generate a Singularity image

To generate an image:

1. Change build args in the [docker-compose file](./singularity/docker-compose.yml) if necessary
1. cd to the singularity directory and run `./build.bash` (takes some time).

---

## Run a slurm job using the image

To run a firedrake script using the singularity image:

1. Copy `sub_fd-singularity.slurm` to the run directory
1. Set your slurm account in the line `#SBATCH --account=[SET_ACCOUNT_HERE]`
1. Set the location of the singularity image in the line `sif_dir=[SET_IMAGE_DIR_HERE]`
1. Change the name of the script if necessary in the line `script="LAPD-like_simplified_CG.py"`
1. Copy your firedrake script to the run directory
1. Submit the job with `sbatch sub_fd-singularity.slurm`