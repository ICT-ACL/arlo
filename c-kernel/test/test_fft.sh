#!/bin/bash

make fft

if [ $? -ne 0 ]; then
	exit 1
fi

echo "Test FFT"
echo "================"
python test_fft.py fft
./fft fft
rm -f ./data/*.dat

echo
echo "Test iFFT"
echo "================"
python test_fft.py ifft
./fft ifft
rm -f ./data/*.dat

rm -f fft