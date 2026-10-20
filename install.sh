#!/bin/sh

base=$(dirname $(realpath "$0"))
pkg=$base/cpp

echo "\nInstalling MarQu C++ lib..."
build_dir="$pkg/build"
[ -d $build_dir ] || mkdir $build_dir
cd $build_dir
cmake ..
make install
ln -sf "$build_dir/compile_commands.json" "$pkg/compile_commands.json"

echo "\nInstalling MarQu python package..."
pkg=$base/python
pip install --user --break-system-packages -e $pkg || pip install --user -e $pkg
