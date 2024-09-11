#!/bin/bash
repo_root=$( cd -- "$(realpath $( dirname -- "${BASH_SOURCE[0]}" ))/../" &> /dev/null && pwd )

cp_from_subdir()
{
    if [ $# -ne 2 ]; then
        echo "cp_from_subdir expected 2 args but got $#"
        exit 2
    fi
    subdir="$1"
    file_or_dir="$2"
    pattern="${repo_root}/$subdir/$file_or_dir"
    files_arr=($pattern)
    if [ ${#files_arr[@]} -eq 0 ]; then
         echo "cp_from_subdir: No file(s)/dir(s) matching pattern [$pattern]"
         exit 1
    fi

    for pth in "${files_arr[@]}"; do
        if [ -e "${pth}" ]; then
            cp -r "${pth}" "${run_dir}"
        else
            echo "No such file/dir [${pth}], skipping..."
        fi
    done
}

# Parse args
if [ $# -eq 2 ] || [ $# -eq 3 ]; then
    run_name="$1"
    script_name="${2%.py}.py"
    if [ $# -eq 3 ]; then
        cfg_suffix="$3"
        cfg_renamed="${script_name%.py}_config.yml"
    else
        cfg_suffix="*"
    fi
    cfg_names="${script_name%.py}_config${cfg_suffix}.yml"
else
    echo "usage $0 [run lbl] [script_name] <config_suffix>"
    exit 2
fi

# Make run directory
run_dir="${repo_root}/runs/${run_name}"
rm -rf "${run_dir}"
mkdir -p "${run_dir}"

# Copy in python script, common dir and config(s)
cp_from_subdir "scripts" "$script_name"
cp_from_subdir "scripts" "common"
cp_from_subdir "config" "$cfg_names"

# If only one config file was copied, remove the suffix
if [ -n "$cfg_renamed" ]; then 
    mv "${run_dir}/${cfg_names}" "${run_dir}/${cfg_renamed}"
fi

echo "Created run dir at $run_dir"
