#!/bin/bash

# script to build pyrad, pyart and PyTDA
# fvj 30.09.2019

if [[ -z "${PYRAD_PATH}" ]]; then
  read -p "Enter the path to your pyrad main directory [default: $HOME/pyrad/]: " PYRAD_PATH
  PYRAD_PATH=${PYRAD_PATH:-$HOME/pyrad/}
  export PYRAD_PATH=$PYRAD_PATH
fi

echo 'Building Pyart...'
./make_pyart.sh

echo 'Building PyTDA...'
./make_pytda.sh

echo 'Building Pyrad...'
./make_pyrad.sh