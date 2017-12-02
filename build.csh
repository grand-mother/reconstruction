#!/bin/csh

echo "Building port_i ..."
cd port_i
make ${1}
echo "--> done"

echo "Building reco ..."
cd ..
make ${1}
echo "--> done"

echo "--> All done"
