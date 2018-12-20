#!/bin/bash

make params

if [ $? -ne 0 ]; then
	exit 1
fi

DATA_DIR='./data'
NROWS=479325
LOW='0.8e8'
HIGH='1.2e8'
LEN=7

if [ ! -d $DATA_DIR ]; then
    mkdir $DATA_DIR
fi

ARGS="--data_dir=$DATA_DIR --nrows=$NROWS --low=$LOW --high=$HIGH --len=$LEN"

echo "Test Params"
echo "================"
for arg in $ARGS; do
    echo $arg
done
echo "----------------"
python test_params.py $ARGS
./params $ARGS
rm -rf $DATA_DIR

rm -f params