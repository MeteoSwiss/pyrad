#!/bin/bash -l
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --ntasks-per-node=1
#SBATCH --mem=64g
#SBATCH --partition=postproc
#SBATCH --account=msrad

# other sbatch options:
# mem-per-cpu : maximum amount of real memory per allocated cpu required by the job
#               --mem >= --mem-per-cpu if --mem is specified

EXECUTABLE=$1

echo "Launched executable $EXECUTABLE"

source activate pyrad

$EXECUTABLE


source deactivate
