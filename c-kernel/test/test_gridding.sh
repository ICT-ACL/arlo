#!/bin/bash

make gridding

if [ $? -ne 0 ]; then
	exit 1
fi

DATA_DIR=./data
VSHAPE=479325,1
GSHAPE=7,1,1024,1024
KSHAPE=1,8,8,8,8

if [ ! -d $DATA_DIR ]; then
    mkdir $DATA_DIR
fi

ARGS="--data_dir=$DATA_DIR --vshape=$VSHAPE --gshape=$GSHAPE --kshape=$KSHAPE"

echo "Test Degridding"
echo "================"
for arg in $ARGS; do
    echo $arg
done
echo "----------------"
python test_gridding.py --mode=degrid $ARGS
./gridding --mode=degrid $ARGS
rm -rf $DATA_DIR/*

echo
echo "Test Gridding"
echo "================"
for arg in $ARGS; do
    echo $arg
done
echo "----------------"
python test_gridding.py --mode=grid $ARGS
./gridding grid --mode=grid $ARGS
rm -rf $DATA_DIR

rm -f gridding