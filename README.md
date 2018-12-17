# Optimized SKA-ARL

## Introduction
- [Algorithm Reference Library (ARL)](https://github.com/SKA-ScienceDataProcessor/algorithm-reference-library)
is a set of algorithms developed for aperture 
synthesis imaging based on Numpy. ARLO is an optimized version of ARL by 
converting computational intensive kernels from python to c++ and applying 
further optimization.

- Applied Computing Lab., ICT, CAS


## Installation
1. install required python packages:
    ```bash
    pip3 install -r requirements.txt 
    ```
2. compile C++ code to `.so` files which can be import as a module in Python code.
    - the C++ code is in `c-kernel/src`
    - there are 4 module should be compiled:
        - [gridding](./c-kernel/src/fourier_transforms): `convolution_gridding_mod.so` by **Makefile.gridding**
        - [fft](./c-kernel/src/fourier_transforms): `fft_support_mod.so` by **Makefile.fft**
        - [cleaners](./c-kernel/src/image): `cleaners_mod.so` by **Makefile.cleaners**
        - [parameters](./c-kernel/src/imaging): `parans_mod.so` by **Makefile.params**
    - use Makefiles to compile: `make -f Makefile.xxx`    
3. move `.so` files to corresponding directories in `arl-python/arl`, e.g.
   ```bash
   mv c-kernel/src/fourier_transforms/convolution_gridding_mod.so arl-python/arl/fourier_transforms
   ```
   similar for other `.so` files.

   
## Demo
- Two notebook demos can be run and test: [imaging-demo1](arl-python/examples/arl/imaging-demo1.ipynb), [imaging-demo2](arl-python/examples/arl/imaging-demo2.ipynb) (modified from the original imaging notebook).
- The intermediate data used by the notebooks are in `arl-python/example/arl/results`.


## Test
`c-kernel/test` is code to test and verify correctness of the optimized C++ code.  