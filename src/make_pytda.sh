#!/bin/bash

# script to build pytda
# fvj 17.11.2016
# sue 28.03.2018

# remove previous built
echo 'Removing previous built...'

if [[ -z "${PYRAD_PATH}" ]]; then
  read -p "Enter the path to your pyrad main directory [default: $HOME/pyrad/]: " PYRAD_PATH
  PYRAD_PATH=${PYRAD_PATH:-$HOME/pyrad/}
  export PYRAD_PATH=$PYRAD_PATH
fi

cd $PYRAD_PATH/src/PyTDA
python setup.py clean --all 

# recompile
echo 'compiling...'
cd $HOME/pyrad/src/PyTDA
python setup.py install