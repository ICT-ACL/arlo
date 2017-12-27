# test/make.inc
# Author: You Haihang, Yang Runkai, Liu Tao,
# Applied Computing Lab, ICT, CAS 

BLASDIR=$(SOFTWARE_DIR)/cblas-3.8.0

INCLUDE_DIR=../src
LIB_DIRS=-L../lib \
		-L$(BLASDIR)/lib

INCLUDES=-I$(INCLUDE_DIR)/utils \
         -I$(INCLUDE_DIR)/fourier_transforms \
         -I$(INCLUDE_DIR)/imaging \
         -I$(INCLUDE_DIR)/image \
         -I$(BLASDIR)/include

CC=g++
CFLAGS=-O3 $(INCLUDES) -fopenmp -DDEBUG
LDFLAGS=$(LIB_DIRS)
LDLIBS=../lib/libfftw3.a $(BLASDIR)/lib/cblas_LINUX.a $(BLASDIR)/lib/blas_LINUX.a -lm -lrt
