#!/bin/sh

pkg=$(dirname $(realpath "$0"))

echo "\nUninstall QuantumCTMC..."
build_dir="$pkg/build"
[ -d $build_dir ] || mkdir $build_dir
cd $build_dir
cmake .. 
make uninstall
