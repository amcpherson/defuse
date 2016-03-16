#!/bin/bash

# Build tools
export CPLUS_INCLUDE_PATH=$PREFIX/src/boost/
cd tools
make
cd ..

# Move to opt directory
mkdir -p $PREFIX/opt/defuse/
cp -R * $PREFIX/opt/defuse/
