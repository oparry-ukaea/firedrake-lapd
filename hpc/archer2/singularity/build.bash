#!/bin/bash

check_ret() {
    ret_code="$?"
    if [ $# -ne 1 ]; then
        echo "check_ret called with $# args (rather than 1) "
        exit 99 
    fi
    lbl=$1
    if [ $ret_code -ne 0 ]; then
        echo "$lbl failed with return code $ret_code"
        exit $ret_code
    fi
}

for cmd in "docker" "singularity"
do
    if ! command -v "$cmd" &> /dev/null
    then
        echo "$cmd command not on path or doesn't exist"
        exit 1
    fi
done

echo "N.B. Pass --no-cache to this script to disable layer caching and force all images to rebuild from scratch"

# Build the docker images
compose_cmd="docker compose build $*"
eval "$compose_cmd"
check_ret "docker compose build"

# Convert to sandbox singularity
singularity build --force \
   --sandbox ./firedrake-singularity \
   docker-daemon://firedrake:latest
check_ret "singularity sandbox build"

# Build singularity image
singularity build --force \
firedrake.sif \
./firedrake-singularity
check_ret "singularity image build"