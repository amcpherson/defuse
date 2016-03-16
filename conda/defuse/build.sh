#!/bin/bash

export CPLUS_INCLUDE_PATH=$PREFIX/src/boost/
cd tools
make
mkdir -p $PREFIX/bin
cp clustermatepairs setcover calccov estislands dosplitalign evalsplitalign localalign splitseq matealign bamfastq $PREFIX/bin


