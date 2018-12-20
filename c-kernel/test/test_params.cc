// test/test_params.cc
// tests of parameters
// Author: You Haihang, Yang Runkai, Liu Tao,
// Applied Computing Lab, ICT, CAS

#include <iostream>
#include <string>
#include <ctime>
#include <getopt.h>

#include "params.h"
#include "utils.h"

using namespace std;

void test_params(string &data_dir, int nrows, int n_unique) {
    double *frequency = new double[nrows];
    double *ufrequency = new double[n_unique];
    int *vmap = new int[nrows];
    int *vmap_exact = new int[nrows];

    load_data(data_dir + "/" + "frequency.dat", frequency, nrows);
    load_data(data_dir + "/" + "ufrequency.dat", ufrequency, n_unique);
    load_data(data_dir + "/" + "vmap.dat", vmap_exact, nrows);

    clock_t start = clock();
    get_rowmap(vmap, frequency, ufrequency, nrows, n_unique);
    clock_t stop = clock();
    printf("Optimized Time: %.2lfs\n", 1.0 * (stop - start) / CLOCKS_PER_SEC);

    double err = diff(vmap, vmap_exact, nrows);
    cout << "Params Error: " << err << endl;

    delete [] frequency;
    delete [] ufrequency;
    delete [] vmap_exact;
    delete [] vmap;
}

int main(int argc, char *argv[]) {
    string data_dir = "./data";
    int nrows = 0;
    double low = 0.0;
    double high = 0.0;
    int len = 0;

    int c;
    while (true) {
        int option_index = 0;
        static struct option long_options[] = {
            {"data_dir", required_argument, 0, 1},
            {"nrows", required_argument, 0, 2},
            {"low", required_argument, 0, 3},
            {"high", required_argument, 0, 4},
            {"len", required_argument, 0, 5}
        };

        c = getopt_long(argc, argv, "", long_options, &option_index);
        if (c == -1) break;

        switch (c) {
            case 1: data_dir = string(optarg); break;
            case 2: nrows = atoi(optarg); break;
            case 3: low = atof(optarg); break;
            case 4: high = atof(optarg); break;
            case 5: len = atoi(optarg); break;
            default: cout << "Wrong options." << endl; break;
        }
    }

    test_params(data_dir, nrows, len);
}
