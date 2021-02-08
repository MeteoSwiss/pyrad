#!/bin/bash

# script to generate the pyrad library reference and remove intermediate outputs
# fvj 20.10.2016

# if command to make pdf from latex exists create latex
latexmk_exists=$(which latexmk 2>/dev/null)
if [ -n "$latexmk_exists" ]; then

    # library reference for developers
    rm -f ../pyrad_library_reference_dev.pdf
    cd library_reference_developers
    make clean
    make latexpdf
    cp build/latex/*.pdf ../../.
    rm -rf build/*
    cd ..

    # library reference for users
    rm -f ../pyrad_library_reference_users.pdf
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
    rm -fr ../../../docs
    git rm -fr ../../../docs
    mkdir ../../../docs
    mv build/html/* ../../../docs
    rm -rf build/*
    touch ../../../docs/.nojekyll
    cd ..
fi
