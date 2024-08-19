#!/bin/bash
repo_root=$( cd -- "$(realpath $( dirname -- "${BASH_SOURCE[0]}" ))/../../" &> /dev/null && pwd )

if [ $# -ne 2 ]; then
    echo "usage $0 [run lbl] [script_name]"
    exit 1
fi
run_name="$1"
script_name="${2%.py}.py"

if [ -z "$ACCOUNT" ]; then
    echo "run export ACCOUNT=your_archer2_account_name before using this script"
    exit 2
fi

script_plus_dependencies="$script_name common"

# Settings
run_dir="${repo_root}/runs/${run_name}"
archer_dir="${repo_root}/hpc/archer2"
rm -rf "${run_dir}"
mkdir -p "${run_dir}"
# Copy Slurm script, setting account, image directory
sed "${archer_dir}/sub_fd-singularity.slurm" -e "s|\[SET_ACCOUNT_HERE\]|$ACCOUNT|" -e "s|\[SET_IMAGE_DIR_HERE\]|${repo_root}/hpc/archer2/singularity|"  -e "s|\[SET_SCRIPT_NAME_HERE\]|${script_name}|"> "${run_dir}/sub_fd-singularity.slurm"
# Copy python scripts and dependencies
for file_or_dir in $script_plus_dependencies; do
    pth="${repo_root}/scripts/$file_or_dir"
    if [ -e "${pth}" ]; then
        cp -r "${pth}" "${run_dir}"
    else
        echo "No such file [${pth}], skipping..."
    fi
done

echo Done
