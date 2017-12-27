// utils/fftw.h
// fftw utils
// Author: You Haihang, Yang Runkai, Liu Tao,
// Applied Computing Lab, ICT, CAS 
#ifndef FFTW_H
#define FFTW_H

#include <complex.h>
#include <fftw3.h>

typedef double complex complex_t;

#define OPTIMIZE_FFTW 0

void fftw_fft2(complex_t *out, complex_t* in, const int m, const int n);

void fftw_ifft2(complex_t *out, complex_t* in, const int m, const int n);

void fftw_fft2(complex_t *out, complex_t *in, const int m, const int n, int backwards);

void fftshift(complex_t *out, complex_t* in, const int m, const int n);

void ifftshift(complex_t *out, complex_t* in, const int m, const int n);

#endif
