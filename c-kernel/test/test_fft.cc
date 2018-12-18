// test/test_fft.cc
// tests of fft
// Author: You Haihang, Yang Runkai, Liu Tao,
// Applied Computing Lab, ICT, CAS

#include <iostream>
#include <string>
#include <ctime>

#include "fft_support.h"
#include "utils.h"

using namespace std;

void test_fft(string &data_dir) {
    int shape[] = {7, 2, 1024, 1024};
    int len = shape[0] * shape[1] * shape[2] * shape[3];
    complex_t *a = new complex_t[len];
    complex_t *af = new complex_t[len];
    complex_t *af_exact = new complex_t[len];

    load_data(data_dir + "/" + "a.dat", a, len);
    load_data(data_dir + "/" + "ia.dat", af_exact, len);

    clock_t start = clock();
    fft(af, a, shape, 4);
    clock_t stop = clock();
    double wtime = 1.0 * (stop - start) / CLOCKS_PER_SEC;
    printf("Optimized Time: %.2lfs\n", wtime);

    double err = diff(af, af_exact, len);
    cout << "FFT Error: " << err << endl;

    delete [] a;
    delete [] af;
    delete [] af_exact;
}


void test_ifft(string &data_dir) {
    int shape[] = {7, 2, 1024, 1024};
    int len = shape[0] * shape[1] * shape[2] * shape[3];
    complex_t *af = new complex_t[len];
    complex_t *a = new complex_t[len];
    complex_t *a_exact = new complex_t[len];

    load_data(data_dir + "/" + "ia.dat", af, len);
    load_data(data_dir + "/" + "a.dat", a_exact, len);

    clock_t start = clock();
    ifft(a, af, shape, 4);
    clock_t stop = clock();
    double wtime = 1.0 * (stop - start) / CLOCKS_PER_SEC;
    printf("Optimized Time: %.2lfs\n", wtime);

    double err = diff(a, a_exact, len);
    cout << "iFFT Error: " << err << endl;

    delete [] af;
    delete [] a;
    delete [] a_exact;
}


int main(int argc, char *argv[]) {
    string data_dir = "./data";

    if (!strcmp(argv[1], "fft")) {
        test_fft(data_dir);
    }
    else if (!strcmp(argv[1], "ifft")) {
        test_ifft(data_dir);
    }

    return 0;
}