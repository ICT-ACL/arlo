#!/bin/bash

make params

if [ $? -ne 0 ]; then
	exit 1
fi

echo "Test Params"
echo "================"
python test_params.py
./params
rm -f ./data/*.dat

rm -f params