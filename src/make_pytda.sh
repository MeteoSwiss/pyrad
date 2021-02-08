#!/bin/bash

# script to build pytda
# fvj 17.11.2016
# sue 28.03.2018

# remove previous built
echo 'Removing previous built...'

cd $HOME/pyrad/src/PyTDA
python setup.py clean --all 

# recompile
echo 'compiling...'
cd $HOME/pyrad/src/PyTDA
python setup.py install