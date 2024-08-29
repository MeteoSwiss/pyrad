#!/bin/bash

# script to build the MeteoSwiss Pyart
# fvj 17.11.2016
# sue 28.03.2018

# fetch python version in use
py=$(python --version)
pyvers=${py:7:3}

# remove previous built
echo 'Removing previous built...'
if [[ -z "${PYRAD_PATH}" ]]; then
  read -p "Enter the path to your pyrad main directory [default: $HOME/pyrad/]: " PYRAD_PATH
  PYRAD_PATH=${PYRAD_PATH:-$HOME/pyrad/}
  export PYRAD_PATH=$PYRAD_PATH
fi

rm -r $HOME/.local/lib/python${pyvers}/site-packages/pyart
rm $HOME/.local/lib/python${pyvers}/site-packages/arm_pyart-*

# recompile
echo 'compiling'
cd $PYRAD_PATH/src/pyart
python setup.py install --user