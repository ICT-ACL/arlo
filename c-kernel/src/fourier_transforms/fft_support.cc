// fourier_transforms/fft_support.cc
// fft functions
// Author: You Haihang, Yang Runkai, Liu Tao,
// Applied Computing Lab, ICT, CAS

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "fft_support.h"

/**
 * Fourier transformation from image to grid space
 * .. note::
 *     If there are four axes then the last outer axes are not transformed
 * :param a: image in `lm` coordinate space
 * :return: `uv` grid
 */
void fft(complex_t *af, complex_t* a, const int *shape, const int n_shapes) {
    assert(n_shapes == 4 || n_shapes == 2);

    if (n_shapes == 4) {
        int m = shape[0], n = shape[1], mm = shape[2], nn = shape[3];
        int len = mm * nn;
        // printf("%d, %d, %d, %d: %d\n", m, n, mm, nn, len);

        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                complex_t *local = (complex_t*)malloc(sizeof(complex_t) * len);
                ifftshift(local, a + (i * n + j) * len, mm, nn);
                fftw_fft2(local, local, mm, nn);
                fftshift(af + (i * n + j) * len, local, mm, nn);
                free(local);
            }
        }
    }

    else {
        int m = shape[0], n = shape[1];
        complex_t *tmp = (complex_t*)malloc(sizeof(complex_t) * m * n);
        ifftshift(af, a, m, n);
        fftw_fft2(tmp, af, m, n);
        fftshift(af, tmp, m, n);
        free(tmp);
    }
}


/**
 * Fourier transformation from image to grid space
 * .. note::
 *     If there are four axes then the last outer axes are not transformed
 * :param a: image in `lm` coordinate space
 * :return: `uv` grid
 */
void ifft(complex_t *a, complex_t* af, const int *shape, const int n_shapes) {
    assert(n_shapes == 4 || n_shapes == 2);

    if (n_shapes == 4) {
        int m = shape[0], n = shape[1], mm = shape[2], nn = shape[3];
        int len = mm * nn;

        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                complex_t *local = (complex_t*)malloc(sizeof(complex_t) * len);
                ifftshift(local, af + (i * n + j) * len, mm, nn);
                fftw_ifft2(local, local, mm, nn);
                fftshift(a + (i * n + j) * len, local, mm, nn);
                free(local);
            }
        }
    }

    else {
        int m = shape[0], n = shape[1];
        complex_t *tmp = (complex_t*)malloc(sizeof(complex_t) * m * n);
        ifftshift(a, af, m, n);
        fftw_ifft2(tmp, a, m, n);
        fftshift(a, tmp, m, n);

        for (int i = 0; i < m * n; i++) {
            a[i] /= (m * n);
        }

        free(tmp);
    }
}


/**
 * Pad a far field image with zeroes to make it the given size.
 * Effectively as if we were multiplying with a box function of the
 * original field's size, which is equivalent to a convolution with a
 * sinc pattern in the uv-grid.
 * .. note::
 *     If there are four axes then the last outer axes are not transformed
 * :param ff: The input far field. Should be smaller than NxN.
 * :param npixel:  The desired far field size
 */
void pad_mid(complex_t *out, complex_t *in, const int npixel, const int *shape, const int n_shapes) {
    assert(n_shapes == 2 || n_shapes == 4);

    int nchan = 1, npol = 1;
    int ny, nx;
    if (n_shapes == 4) {
        nchan = shape[0], npol = shape[1], ny = shape[2], nx = shape[3];
    }
    else {
        ny = shape[0], nx = shape[1];
    }

    int n = nchan * npol * npixel * npixel;

    if (npixel == nx) {
        memcpy(out, in, sizeof(complex_t) * n);
    }

    assert(nx == ny && npixel > nx);

    int sx = (npixel - nx) / 2, sy = (npixel - ny) / 2;
    memset(out, 0, sizeof(complex_t) * n);

    for (int c = 0; c < nchan; c++) {
        for (int p = 0; p < npol; p++) {
            for (int i = 0; i < ny; i++) {
                memcpy(out + (((c * npol) + p) * npixel + (sy + i)) * npixel + sx,
                       in + (((c * npol) + p) * ny + i) * nx,
                       sizeof(complex_t) * nx);
            }
        }
    }

}


/**
 * Extract a section from middle of a map
 * Suitable for zero frequencies at npixel/2. This is the reverse
 * operation to pad.
 * .. note::
 *     If there are four axes then the last outer axes are not transformed
 * :param npixel:
 * :param a: grid from which to extract
 */
void extract_mid(complex_t *out, complex_t *in, const int npixel, const int *shape, const int n_shapes) {
    assert(n_shapes == 2 || n_shapes == 4);

    int nchan = 1, npol = 1;
    int ny, nx;
    if (n_shapes == 4) {
        nchan = shape[0], npol = shape[1], ny = shape[2], nx = shape[3];
    }
    else {
        ny = shape[0], nx = shape[1];
    }

    int n = nchan * npol * npixel * npixel;

    if (npixel == nx) {
        memcpy(out, in, sizeof(complex_t) * n);
    }

    assert(nx == ny && npixel < nx);

    int sx = (nx - npixel) / 2, sy = (ny - npixel) / 2;
    memset(out, 0, sizeof(complex_t) * n);

    for (int c = 0; c < nchan; c++) {
        printf("c = %d\n", c);
        for (int p = 0; p < npol; p++) {
            for (int i = 0; i < npixel; i++) {
                memcpy(out + (((c * npol) + p) * npixel + i) * npixel,
                    in + (((c * npol) + p) * ny + (sy + i)) * nx + sx,
                    sizeof(complex_t) * npixel);
            }
        }
    }

}


/**
 * Extract the (xf-th,yf-th) w-kernel from the oversampled parent
 * Offsets are suitable for correcting of fractional coordinates,
 * e.g. an offset of (xf,yf) results in the kernel for an (-xf,-yf)
 * sub-grid offset.
 * We do not want to make assumptions about the source grid's symmetry
 * here, which means that the grid's side length must be at least
 * kernel_oversampling*(npixel+2) to contain enough information in all circumstances
 * :param xf:
 * :param yf:
 * :param a: grid from which to extract
 * :param kernel_oversampling: oversampling factor
 * :param kernelwidth: size of section
 */
void extract_oversampled(complex_t *mid, complex_t *a, int xf, int yf, int kernel_oversampling, int kernelwidth, const int *shape) {
    assert(xf >= 0 && xf < kernel_oversampling);
    assert(yf >= 0 && yf < kernel_oversampling);

    int p = shape[0], q = shape[1];
    int npixela = p;
    int my = npixela / 2 - kernel_oversampling * (kernelwidth / 2) - yf;
    int mx = npixela / 2 - kernel_oversampling * (kernelwidth / 2) - xf;

    assert(mx >= 0 && my >= 0);

    for (int i = 0; i < kernelwidth; i++) {
        for (int j = 0; j < kernelwidth; j++) {
            mid[i*kernelwidth+j] = kernel_oversampling * kernel_oversampling *
                a[(my + kernel_oversampling * i) * q + (mx + kernel_oversampling * j)];
        }
    }
}



int fft_c(PyArrayObject *&af_obj, PyArrayObject *&a_obj) {

    const int n_shapes = a_obj->nd;
    int shape[n_shapes];
    for (int i = 0; i < n_shapes; i++) {
        shape[i] = a_obj->dimensions[i];
    }

    complex_t *af = (complex_t*)af_obj->data;
    complex_t *a = (complex_t*)a_obj->data;

    fft(af, a, shape, n_shapes);

    return 0;
}


int ifft_c(PyArrayObject *&a_obj, PyArrayObject *&af_obj) {

    const int n_shapes = af_obj->nd;
    int shape[n_shapes];
    for (int i = 0; i < n_shapes; i++) {
        shape[i] = af_obj->dimensions[i];
    }

    complex_t *a = (complex_t*)a_obj->data;
    complex_t *af = (complex_t*)af_obj->data;

    ifft(a, af, shape, n_shapes);

    return 0;
}


int pad_mid_c(PyArrayObject *&out_obj, PyArrayObject *&in_obj, int npixel) {

    const int n_shapes = in_obj->nd;

    int shape[n_shapes];
    for (int i = 0; i < n_shapes; i++) {
        shape[i] = in_obj->dimensions[i];
    }

    complex_t *out = (complex_t*)out_obj->data;
    complex_t *in = (complex_t*)in_obj->data;

    pad_mid(out, in, npixel, shape, n_shapes);

    return 0;
}


int extract_mid_c(PyArrayObject *&out_obj, PyArrayObject *&in_obj, int npixel) {

    const int n_shapes = in_obj->nd;

    int shape[n_shapes];
    for (int i = 0; i < n_shapes; i++) {
        shape[i] = in_obj->dimensions[i];
    }

    complex_t *out = (complex_t*)out_obj->data;
    complex_t *in = (complex_t*)in_obj->data;

    pad_mid(out, in, npixel, shape, n_shapes);

    return 0;
}


int extract_oversampled_c(PyArrayObject *&mid_obj, PyArrayObject *&a_obj,
        int xf, int yf, int kernel_oversampling, int kernelwidth) {

    const int n_shapes = a_obj->nd;

    int shape[n_shapes];
    for (int i = 0; i < n_shapes; i++) {
        shape[i] = a_obj->dimensions[i];
    }

    complex_t *mid = (complex_t*)mid_obj->data;
    complex_t *a = (complex_t*)a_obj->data;

    extract_oversampled(mid, a, xf, yf, kernel_oversampling, kernelwidth, shape);

    return 0;
}
