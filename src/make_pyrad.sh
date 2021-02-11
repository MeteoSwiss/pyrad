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

# clean pyrad
echo 'icleaning build..'
cd $HOME/pyrad/src/pyrad_proc
python setup.py clean --all

# recompile
echo 'compiling...'
cd $HOME/pyrad/src/pyrad_proc
python setup.py install --user