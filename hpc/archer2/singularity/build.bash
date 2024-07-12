#!/bin/bash

# Check docker and singularity are available
for cmd in "docker" "singularity"
do
    if ! command -v "$cmd" &> /dev/null
    then
        echo "$cmd command not on path or doesn't exist"
        exit 1
    fi
done

# Build the docker images
docker compose build

# Convert to sandbox singularity
singularity build --force \
   --sandbox ./firedrake-singularity \
   docker-daemon://firedrake:latest

# Build singularity image
singularity build --force \
firedrake.sif \
./firedrake-singularity