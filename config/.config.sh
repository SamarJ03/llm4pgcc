#!/bin/bash

conda env remove -n .pgccEnv

env_name=".pgccEnv"

eval "$(conda shell.bash hook)"
conda env create -p "./$env_name" -f config/environment.yml python=3.10.13
# conda rename -n "$env_name" .pgccEnv

conda activate "./$env_name"
# source activate pgccEnv