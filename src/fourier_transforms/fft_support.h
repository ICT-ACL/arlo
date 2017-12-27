#ifndef FFT_SUPPORT_H
#define FFT_SUPPORT_H

#include "../utils/fftw.h"


void fft(complex_t *af, complex_t* a, const int *shape, const int n_shapes);

void ifft(complex_t *a, complex_t* af, const int *shape, const int n_shapes);

void pad_mid(complex_t *out, complex_t *in, const int npixel, const int *shape, const int n_shapes);

void extract_mid(complex_t *out, complex_t *in, const int npixel, const int *shape, const int n_shapes);

void extract_oversampled(complex_t *mid, complex_t *a, int xf, int yf, int kernel_oversampling, int kernelwidth, const int *shape);

#endif