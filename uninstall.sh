#!/bin/sh

base=$(dirname $(realpath "$0"))
pkg=$base/cpp

echo "\nUninstalling MarQu C++ lib..."
build_dir="$pkg/build"
[ -d $build_dir ] || mkdir $build_dir
cd $build_dir
cmake .. 
make uninstall

echo "\nUninstalling MarQu python package..."
pkg=$base/python
pip uninstall --break-system-packages marqu || pip uninstall marqu
