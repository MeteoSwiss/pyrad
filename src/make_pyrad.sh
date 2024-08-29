#!/bin/bash

# script to build pyrad
# fvj 17.11.2016
# sue 28.03.2018

# fetch python version in use
py=$(python --version)
pyvers=${py:7:3}

# remove previous built
echo 'Removing previous built...'
 
rm -r $HOME/.local/lib/python${pyvers}/site-packages/pyrad
rm -r $HOME/.local/lib/python${pyvers}/site-packages/mch_pyrad-*

if [[ -z "${PYRAD_PATH}" ]]; then
  read -p "Enter the path to your pyrad main directory [default: $HOME/pyrad/]: " PYRAD_PATH
  PYRAD_PATH=${PYRAD_PATH:-$HOME/pyrad/}
    export PYRAD_PATH=$PYRAD_PATH
fi

# clean pyrad
echo 'icleaning build..'
cd $PYRAD_PATH/src/pyrad_proc
python setup.py clean --all

# recompile
echo 'compiling...'
cd $PYRAD_PATH/src/pyrad_proc
python setup.py install --user
