#!/bin/bash

# Build tools
export CPLUS_INCLUDE_PATH=$PREFIX/src/boost/
cd tools
make
cd ..

# Move to opt directory
mkdir -p $PREFIX/opt/defuse/
cp -R * $PREFIX/opt/defuse/

# Create symlinks for defuse and create_reference_dataset
ln -s $PREFIX/opt/defuse/scripts/defuse_run.pl $PREFIX/bin/defuse_run.pl
ln -s $PREFIX/opt/defuse/scripts/defuse_create_ref.pl $PREFIX/bin/defuse_create_ref.pl
ln -s $PREFIX/opt/defuse/scripts/defuse_get_reads.pl $PREFIX/bin/defuse_get_reads.pl
