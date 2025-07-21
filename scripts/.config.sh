#!/bin/bash

env_name=".pgccEnv"
cwd=""
if [ -n "$1" ]; then
    cwd="$1"

eval "$(conda shell.bash hook)"
conda env create -p "./$env_name" -f config/environment.yml python=3.10.13
conda rename -p "$cwd/$env_name" .pgccEnv

conda activate $env_name
# source activate pgccEnv