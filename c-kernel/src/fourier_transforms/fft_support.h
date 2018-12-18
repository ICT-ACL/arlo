// fourier_transforms/fft_support.h
// fft functions
// Author: You Haihang, Yang Runkai, Liu Tao,
// Applied Computing Lab, ICT, CAS 

#ifndef FFT_SUPPORT_H
#define FFT_SUPPORT_H

#include <ndarrayobject.h>
#include "fftw.h"

void fft(complex_t *af, complex_t* a, const int *shape, const int n_shapes);
void ifft(complex_t *a, complex_t* af, const int *shape, const int n_shapes);
void pad_mid(complex_t *out, complex_t *in, const int npixel, const int *shape, const int n_shapes);
void extract_mid(complex_t *out, complex_t *in, const int npixel, const int *shape, const int n_shapes);
void extract_oversampled(complex_t *mid, complex_t *a, int xf, int yf, int kernel_oversampling, int kernelwidth, const int *shape);

int fft_c(PyArrayObject *&af_obj, PyArrayObject *&a_obj);
int ifft_c(PyArrayObject *&a_obj, PyArrayObject *&af_obj);
int pad_mid_c(PyArrayObject *&out_obj, PyArrayObject *&in_obj, int npixel);
int extract_mid_c(PyArrayObject *&out_obj, PyArrayObject *&in_obj, int npixel);
int extract_oversampled_c(PyArrayObject *&mid_obj, PyArrayObject *&a_obj,
        int xf, int yf, int kernel_oversampling, int kernelwidth);

#endif
