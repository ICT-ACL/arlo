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
1. Compile code in `src` first and move the `.so` libs to corresponding directories to make sure some python modules imported correctly.
2. Run scripts to verify the results of original python code and optimized c++ code.
    - For all tests, run `sh test_all.sh`.
    - For params, run `sh test_params.sh`.
    - For gridding, run `sh test_gridding.sh`.
    - For fft, run `sh test_fft.sh`.  
    - For cleaners, run `sh test_cleaners.sh`.
3. The data in test code is randomly created.
