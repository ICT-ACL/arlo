#!/bin/bash

make fft

if [ $? -ne 0 ]; then
	exit 1
fi

DATA_DIR=./data
SHAPE=7,2,1024,1024

if [ ! -d $DATA_DIR ]; then
    mkdir $DATA_DIR
fi

ARGS="--data_dir=$DATA_DIR --shape=$SHAPE"

echo "Test FFT"
echo "================"
for arg in $ARGS; do
    echo $arg
done
echo "----------------"
python test_fft.py --mode=fft $ARGS
./fft --mode=fft $ARGS
rm -rf $DATA_DIR/*

echo
echo "Test iFFT"
echo "================"
for arg in $ARGS; do
    echo $arg
done
echo "----------------"
python test_fft.py --mode=fft $ARGS
./fft --mode=fft $ARGS
rm -rf $DATA_DIR

rm -f fft