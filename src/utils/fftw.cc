#include <map>
#include <string.h>
#include "fftw.h"                         

std::map<std::pair<int,int>, fftw_plan> fftw_plans;

// fft2d
void fftw_fft2(complex_t *out, complex_t* in, const int m, const int n) {
    fftw_fft2(out, in, m, n, 0);
}


// inverse fft2d
void fftw_ifft2(complex_t *out, complex_t* in, const int m, const int n) {
    fftw_fft2(out, in, m, n, 1);
}


// fft2d function
void fftw_fft2(complex_t *out, complex_t *in, const int m, const int n, int backwards){

    fftw_plan p;

    if (OPTIMIZE_FFTW) {
        std::pair<int,int> shape = std::make_pair(m, n);

        if (fftw_plans.find(shape) == fftw_plans.end()) {

            complex_t *in      = (complex_t *)fftw_malloc(sizeof(*in)      * m * n);
            complex_t *out_tmp = (complex_t *)fftw_malloc(sizeof(*out_tmp) * m * n);

            p = fftw_plan_dft_2d(m, n, in, out_tmp, 
                                 backwards ? FFTW_BACKWARD:FFTW_FORWARD,
                                 FFTW_ESTIMATE);

            fftw_plans.insert(std::make_pair(shape, p));

            fftw_free(in);
            fftw_free(out_tmp);
        }
    }   

    p = fftw_plan_dft_2d(m, n, in, out, 
                         backwards ? FFTW_BACKWARD:FFTW_FORWARD,
                         FFTW_ESTIMATE);

    if (backwards) {
        double coeff = 1.0 / (m * n);
        for (int i = 0; i < m * n; i++) {
            out[i] *= coeff;
        }
    }

    fftw_execute(p);
    fftw_destroy_plan(p);
}


// fft shift
// notice that in and out cannot be the same array
void fftshift(complex_t *out, complex_t* in, const int m, const int n) {
    // 1 --> 4
    for (int i = 0; i < m - m/2; i++) {
        memcpy(out + (m/2 + i) * n + (n/2), in + i * n, sizeof(complex_t) * (n - n/2));
    }

    // 2 --> 3
    for (int i = 0; i < m - m/2; i++) {
        memcpy(out + (m/2 + i) * n, in + i * n + (n - n/2), sizeof(complex_t) * (n/2));
    }

    // 3 --> 2
    for (int i = 0; i < m/2; i++) {
        memcpy(out + i * n + (n/2), in + (i + m - m/2) * n, sizeof(complex_t) * (n - n/2));
    }

    // 4 --> 1
    for (int i = 0; i < m/2; i++) {
        memcpy(out + i * n, in + (i + m - m/2) * n + (n - n/2), sizeof(complex_t) * (n/2));
    }
}


// inverse fft shift
// notice that in and out cannot be the same array
void ifftshift(complex_t *out, complex_t* in, const int m, const int n) {
    // 1 --> 4
    for (int i = 0; i < m/2; i++) {
        memcpy(out + (i + m - m/2) * n + (n - n/2), in + i * n, sizeof(complex_t) * (n/2));
    }

    // 2 --> 3
    for (int i = 0; i < m/2; i++) {
        memcpy(out + (i + m - m/2) * n, in + i * n + (n/2), sizeof(complex_t) * (n - n/2));
    }

    // 3 --> 2
    for (int i = 0; i < m - m/2; i++) {
        memcpy(out + i * n + (n - n/2), in + (m/2 + i) * n, sizeof(complex_t) * (n/2));
    }

    // 4 --> 1
    for (int i = 0; i < m - m/2; i++) {
        memcpy(out + i * n, in + (m/2 + i) * n + (n/2), sizeof(complex_t) * (n - n/2));
    }
}
