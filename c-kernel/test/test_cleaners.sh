#!/bin/bash

make cleaners

if [ $? -ne 0 ]; then
	exit 1
fi

echo "Test Cleaners"
echo "================"
python test_cleaners.py
./cleaners
rm -f ./data/*.dat

rm -f cleaners