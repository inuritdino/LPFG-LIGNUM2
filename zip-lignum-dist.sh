#!/bin/bash

# version of the LPFG-LIGNUM code
VER=2.0.0

# Source files to copy
CPPS="brent.cpp methods.cpp mtg.cpp randist.cpp scatter.cpp shadow_propagation.cpp"
# Header files to copy
HPPS="brent.hpp lgm.h lgmconst.h lgmtypes.h Perttunen1998.par Sievanen2008.par user.h"
# LPFG files
LPFGS="lsystem.l view.v material.mat bark3-128.rgb" 
# Helping files
HELPS="README.md NEWS"

# Create temporary directory
DISTDIR=lignum-dist-$VER
mkdir $DISTDIR

# Copy the file into the directory
cp $CPPS $HPPS $LPFGS $HELPS $DISTDIR

# Create an archive
ZIPARCH=lignum-dist-$VER.zip
zip -r $ZIPARCH $DISTDIR

# Remove the directory
rm -rf $DISTDIR
