#!/bin/bash

# script to generate the MeteoSwiss Pyart library reference and remove intermediate outputs
# fvj 20.10.2016

# if command to make pdf from latex exists create latex
latexmk_exists=$(which latexmk 2>/dev/null)
if [ -n "$latexmk_exists" ]; then

    # library reference for developers
    rm -f ../pyart-mch_library_reference_dev.pdf
    cd library_reference_developers
    make clean
    make latexpdf
    cp build/latex/*.pdf ../../.
    rm -rf build/*
    cd ..

    # library reference for users
    # if command to make pdf from latex exists create latex
    rm -f ../pyart-mch_library_reference_users.pdf
    cd library_reference_users
    make clean
    make latexpdf
    cp build/latex/*.pdf ../../.
    rm -rf build/*
    cd ..
else
    echo 'Unable to produce latex files. Command latexmk not installed'
fi

# if we are in master branch make html for github pages
branch=$(git branch | grep "*" | cut  -d ' ' -f 2)
if [ "${branch}" = "master" ];then
    cd library_reference_users
    make clean
    make html
    cd ../../../src/pyart/
    rm -fr ./docs
    git rm -fr ./docs
    mkdir -p docs
    cd ../../doc/pyart-mch/library_reference_users
    mv build/html/* ../../../src/pyart/docs
    rm -rf build/*
    touch ../../../src/pyart/docs/.nojekyll
    cd ..
fi
