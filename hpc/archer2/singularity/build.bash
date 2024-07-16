#!/bin/bash

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

# Convert to sandbox singularity
singularity build --force \
   --sandbox ./firedrake-singularity \
   docker-daemon://firedrake:latest

# Build singularity image
singularity build --force \
firedrake.sif \
./firedrake-singularity