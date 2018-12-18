# test/Makefile
# Author: You Haihang, Yang Runkai, Liu Tao,
# Applied Computing Lab, ICT, CAS 

.PHONY: clean all

include make.inc

SRC_DIR=../src

SRCS=$(SRC_DIR)/utils/fftw.cc \
	 $(SRC_DIR)/utils/debug.cc \
	 $(SRC_DIR)/fourier_transforms/convolutional_gridding.cc \
	 $(SRC_DIR)/fourier_transforms/fft_support.cc \
	 $(SRC_DIR)/image/cleaners.cc \
	 $(SRC_DIR)/imaging/params.cc

OBJS=$(SRCS:.cc=.o)

%.o: %.cc $(INC)
	$(CC) -c -o $@ $< $(CFLAGS)

all: gridding fft cleaners params

gridding: $(OBJS) test_gridding.cc
	$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS) $(LDLIBS)

fft: $(OBJS) test_fft.cc
	$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS) $(LDLIBS)

cleaners: $(OBJS) test_cleaners.cc
	$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS) $(LDLIBS)

params: $(OBJS) test_params.cc
	$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS) $(LDLIBS)

clean:
	rm -f $(OBJS) gridding fft cleaners params *.o
