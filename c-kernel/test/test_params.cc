// test/test_params.cc
// tests of parameters
// Author: You Haihang, Yang Runkai, Liu Tao,
// Applied Computing Lab, ICT, CAS

#include <iostream>
#include <string>
#include <ctime>

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
    test_params(data_dir, 479325, 7);
}
