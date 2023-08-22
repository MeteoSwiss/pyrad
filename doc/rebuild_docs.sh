# script to rebuild documentation after removing intermediates
rm -r build
rm -r source/API/generated/*
make html
