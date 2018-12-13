# test/Makefile
# Author: You Haihang, Yang Runkai, Liu Tao,
# Applied Computing Lab, ICT, CAS 

.PHONY: clean all

include make.inc

SRC_DIR=../src

SRCS=$(SRC_DIR)/utils/fftw.cc \
	 $(SRC_DIR)/utils/debug.cc \
	 $(SRC_DIR)/fourier_transforms/fft_support.cc \
	 $(SRC_DIR)/fourier_transforms/convolutional_gridding.cc \
	 $(SRC_DIR)/image/cleaners.cc \
	 $(SRC_DIR)/imaging/params.cc 


OBJS=$(SRCS:.cc=.o)
	
TEST=test_clean.cc

%.o: %.cc $(INC) 
	$(CC) -c -o $@ $< $(CFLAGS)

all: test

test: $(OBJS) $(TEST)
	$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS) $(LDLIBS)

clean:
	rm -f $(OBJS) test
