// test/test_fft.cc
// tests of fft
// Author: You Haihang, Yang Runkai, Liu Tao,
// Applied Computing Lab, ICT, CAS

#include <iostream>
#include <string>
#include <vector>
#include <ctime>
#include <getopt.h>

#include "fft_support.h"
#include "utils.h"

using namespace std;


void test_fft(string &data_dir, vector<int> &shape) {
    int len = 1;
    for (int i = 0; i < shape.size(); i++) {
        len *= shape[i];
    }
    complex_t *a = new complex_t[len];
    complex_t *af = new complex_t[len];
    complex_t *af_exact = new complex_t[len];

    load_data(data_dir + "/" + "a.dat", a, len);
    load_data(data_dir + "/" + "ia.dat", af_exact, len);

    clock_t start = clock();
    fft(af, a, &shape[0], shape.size());
    clock_t stop = clock();
    double wtime = 1.0 * (stop - start) / CLOCKS_PER_SEC;
    printf("Optimized Time: %.2lfs\n", wtime);

    double err = diff(af, af_exact, len);
    cout << "FFT Error: " << err << endl;

    delete [] a;
    delete [] af;
    delete [] af_exact;
}


void test_ifft(string &data_dir, vector<int> &shape) {
    int len = 1;
    for (int i = 0; i < shape.size(); i++) {
        len *= shape[i];
    }
    complex_t *af = new complex_t[len];
    complex_t *a = new complex_t[len];
    complex_t *a_exact = new complex_t[len];

    load_data(data_dir + "/" + "ia.dat", af, len);
    load_data(data_dir + "/" + "a.dat", a_exact, len);

    clock_t start = clock();
    ifft(a, af, &shape[0], shape.size());
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
    string mode = "fft";
    string shape_str = "0,0,0,0";

    int c;
    while (true) {
        int option_index = 0;
        static struct option long_options[] = {
            {"data_dir", required_argument, 0, 1},
            {"mode", required_argument, 0, 2},
            {"shape", required_argument, 0, 3}
        };

        c = getopt_long(argc, argv, "", long_options, &option_index);
        if (c == -1) break;

        switch (c) {
            case 1: data_dir = string(optarg); break;
            case 2: mode = string(optarg); break;
            case 3: shape_str = string(optarg); break;
            default: cout << "Wrong options." << endl; break;
        }
    }

    vector<int> shape = str2arr(shape_str);
    if (shape.size() == 2 || shape.size() == 4) {
        if (mode == "fft") {
            test_fft(data_dir, shape);
        }
        else if (mode == "ifft") {
            test_ifft(data_dir, shape);
        }
        else {}
    }

    return 0;
}