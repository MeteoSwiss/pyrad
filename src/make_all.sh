#!/bin/bash

# script to build pyrad, pyart and PyTDA
# fvj 30.09.2019

echo 'Building Pyart...'
./make_pyart.sh

echo 'Building PyTDA...'
./make_pytda.sh

echo 'Building Pyrad...'
./make_pyrad.sh
