#!/bin/bash
#SBATCH -N 1                 # nodes requested
#SBATCH -c 1                 # cores requested
#SBATCH -t 10:0:00           # time requested in hour:minute:second
#SBATCH --partition postproc

# Default values
DEFAULT_PYART_CONFIG="$HOME/pyrad/config/pyart/mch_config.py"
DEFAULT_CONDA_ENV_NAME="pyart"

# Assign arguments or defaults
PYART_CONFIG=${3:-$DEFAULT_PYART_CONFIG}

export PYTHONUNBUFFERED=1
export PATH=$PATH
export PYART_CONFIG=$PYART_CONFIG
export CONDA_ENV_NAME=$DEFAULT_CONDA_ENV_NAME

current_time=$(date +'%Y%m%d%H%M%S')

echo $current_time
echo "Launched executable $EXECUTABLE"

source activate $CONDA_ENV_NAME

$EXECUTABLE

source deactivate



#!/bin/bash
#SBATCH -N 1     # nodes requested
#SBATCH -c 1      # cores requested
#SBATCH -t 10:0:00  # time requested in hour:minute:second
#SBATCH --partition postproc

# other sbatch options:
# mem-per-cpu : maximum amount of real memory per allocated cpu required by the job
#               --mem >= --mem-per-cpu if --mem is specified

EXECUTABLE=$1

export PYTHONUNBUFFERED=1
export PATH=/store_new/mch/msrad/utils/miniconda-wolfensb/envs/pyart/bin/:$PATH
export PYART_CONFIG=$HOME/pyrad/config/pyart/mch_config.py
current_time=$(date +'%Y%m%d%H%M%S')

echo
echo "Launched executable $EXECUTABLE"

source activate pyart

$EXECUTABLE

source deactivate
