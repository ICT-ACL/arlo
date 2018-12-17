#!/bin/bash

make -f Makefile.gridding

if [ $? -ne 0 ]; then
	exit 1
fi

echo "Test Degridding"
echo "================"
python test_gridding.py degrid
./test degrid
rm -f ./data/*.dat

echo

echo "Test Gridding"
echo "================"
python test_gridding.py grid
./test grid
rm -f ./data/*.dat

rm -f test