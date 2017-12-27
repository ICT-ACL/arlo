.PHONY: clean all

PYTHON_DIR=/home/yangrunkai/myinstall/.pyenv/versions/anaconda3-4.4.0
INCLUDE_DIR=../../src
SRC_DIR=../../src
LIB_DIR=../../lib

CC=g++

INCLUDES=-I$(PYTHON_DIR)/include/python3.6m \
		 -I$(PYTHON_DIR)/lib/python3.6/site-packages/numpy/core/include/numpy \
		 -I$(INCLUDE_DIR)/utils \
         -I$(INCLUDE_DIR)/fourier_transforms \
         -I$(INCLUDE_DIR)/imaging \
         -I$(INCLUDE_DIR)/image

LIBS=-L$(PYTHON_DIR)/lib

CFLAGS=-O3 -fopenmp -pthread -fPIC $(INCLUDES)

LDFLAGS=-shared $(LIBS)

LDLIBS=-lpython3.6m -lm -lrt
# $(LIB_DIR)/cblas_LINUX.a $(LIB_DIR)/blas_LINUX.a \


SRCS=$(SRC_DIR)/image/cleaners.cc \
	 $(SRC_DIR)/image/cleaners_mod.cc

OBJS=$(SRCS:.cc=.o)
TARGET=cleaners_mod.so

%.o: %.cc
	$(CC) -c $< -o $@ $(CFLAGS)


all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS) $(LDLIBS)


clean:
	rm -f *.o *.so
