#!/bin/bash

set -e
set -u

source_dir=$PWD/acts
install_dir=$PWD/acts-install
build_dir=$PWD/acts-build

cmake -S $source_dir -B $build_dir -GNinja \
-DCMAKE_BUILD_TYPE=Release  \
-DCMAKE_CXX_COMPILER_LAUNCHER=ccache \
-DCMAKE_INSTALL_PREFIX=$install_dir \
-DACTS_BUILD_EXAMPLES=ON \
-DACTS_BUILD_EXAMPLES_DD4HEP=ON \
-DACTS_BUILD_EXAMPLES_GEANT4=ON \
-DACTS_BUILD_PLUGIN_DD4HEP=ON \
-DACTS_BUILD_ANALYSIS_APPS=ON \
-DACTS_BUILD_EXAMPLES_PYTHON_BINDINGS=ON

cmake --build $build_dir

cmake --install $build_dir
