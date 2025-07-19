#!/bin/sh
# point at your FairSoft v16 install:
source /afs/cern.ch/work/p/pasenov/fairsoft/install/bin/FairRootConfig.sh
export FAIRROOT=/afs/cern.ch/work/p/pasenov/fairsoft/install
export ROOTSYS=$FAIRROOT
export PATH=$ROOTSYS/bin:$PATH
export LD_LIBRARY_PATH=$ROOTSYS/lib64:$LD_LIBRARY_PATH

# for sanity
echo "==> Using ROOT from:    $(which root-config)"
echo "    ROOTSYS =           $ROOTSYS"
echo "    resource-dir =      $(root-config --resource-dir)"

# get ROOT include & libs
ROOTINCDIR=$(root-config --incdir)
ROOTLIBS=$(root-config --cflags --libs)

# FairRoot headers & libs
FAIRROOTINCDIR="$FAIRROOT/include"
FAIRROOTLIBDIR="$FAIRROOT/lib64"

echo "FAIRROOTINCDIR = ${FAIRROOTINCDIR}"
echo "FAIRROOTLIBDIR = ${FAIRROOTLIBDIR}"
echo "ROOTINCDIR     = ${ROOTINCDIR}"
echo "ROOTLIBS       = ${ROOTLIBS}"

OPTIONS="-O3 -Wall"
echo "OPTIONS = ${OPTIONS}"

# compile your tools
g++ ${OPTIONS} order_input.cc  -o order_input    -std=c++17
g++ ${OPTIONS} split_input.cc  -o split_input    -std=c++17

g++ ${OPTIONS} skim.cc Station_hits.h releff.h \
    -o skim \
    -I${ROOTINCDIR} -I${FAIRROOTINCDIR} \
    ${ROOTLIBS} -lX11 -std=c++17

g++ ${OPTIONS} add_output.cc releff.h \
    -o add_output \
    -I${ROOTINCDIR} \
    ${ROOTLIBS} -lX11 -std=c++17
