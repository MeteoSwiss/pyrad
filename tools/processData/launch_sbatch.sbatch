#!/bin/bash
#SBATCH -N 1                 # nodes requested
#SBATCH -c 1                 # cores requested
#SBATCH -t 10:0:00           # time requested in hour:minute:second
#SBATCH --partition postproc

# Default values
DEFAULT_PYART_CONFIG="$HOME/pyrad/config/pyart/mch_config.py"
DEFAULT_CONDA_ENV_NAME="pyart"

# Assign arguments or defaults
EXECUTABLE=$1
PYART_CONFIG=${2:-$DEFAULT_PYART_CONFIG}

export PATH=$PATH
export PYART_CONFIG=$PYART_CONFIG
export CONDA_ENV_NAME=$DEFAULT_CONDA_ENV_NAME

current_time=$(date +'%Y%m%d%H%M%S')

echo $current_time
echo "Launched executable $EXECUTABLE"

$EXECUTABLE

