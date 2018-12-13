# imaging/Makefile
# Makefile of params
# Author: You Haihang, Yang Runkai, Liu Tao,
# Applied Computing Lab, ICT, CAS 
.PHONY: clean all

PYTHON_DIR=/home/yangrunkai/myinstall/.pyenv/versions/anaconda3-4.4.0

CC=g++

INCLUDES=-I$(PYTHON_DIR)/include/python3.6m \
		 -I$(PYTHON_DIR)/lib/python3.6/site-packages/numpy/core/include/numpy
LIBS=-L$(PYTHON_DIR)/lib

CFLAGS=-O3 -pthread -fPIC $(INCLUDES)
LDFLAGS=-shared $(LIBS)
LDLIBS=-lpython3.6m

SRCS=params.cc params_mod.cc
OBJS=$(SRCS:.cc=.o)
TARGET=params_mod.so

%.o: %.cc
	$(CC) -c $< -o $@ $(CFLAGS)


all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS) $(LDLIBS)


clean:
	rm -f *.o *.so
