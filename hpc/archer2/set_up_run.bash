#!/bin/bash
repo_root=$( cd -- "$(realpath $( dirname -- "${BASH_SOURCE[0]}" ))/../../" &> /dev/null && pwd )

if [ $# -ne 1 ]; then
    echo "usage $0 [run lbl"]
    exit 1
fi
run_name="$1"
if [ -z "$ACCOUNT" ]; then
    echo "run export ACCOUNT=your_archer2_account_name before using this script"
    exit 2
fi

script_name="LAPD-like_simplified_CG.py"

# Settings
run_dir="${repo_root}/runs/${run_name}"
archer_dir="${repo_root}/hpc/archer2"
script="${repo_root}/scripts/${script_name}"
rm -rf "${run_dir}"
mkdir -p "${run_dir}"
# Copy Slurm script, setting account, image directory
sed "${archer_dir}/sub_fd-singularity.slurm" -e "s|\[SET_ACCOUNT_HERE\]|$ACCOUNT|" -e "s|\[SET_IMAGE_DIR_HERE\]|${repo_root}/hpc/archer2/singularity|" > "${run_dir}/sub_fd-singularity.slurm"
# Copy python script
cp "$script" "${run_dir}"
