// test/test_gridding.cc
// tests of convolutional gridding
// Author: You Haihang, Yang Runkai, Liu Tao,
// Applied Computing Lab, ICT, CAS

#ifndef UTILS_H
#define UTILS_H

#include <complex.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <algorithm>

typedef double complex complex_t;

using namespace std;

inline void load_line(string &line, int *arr, int len) {
    stringstream sin(line);
    string item;
    int val;

    sin >> item >> item;
    for (int i = 0; i < len; i++) {
        sin >> arr[i];
    }
}

// vshape: shape of vis
// gshape: shape of uvgrid
// kshape: shape of kernels
void load_shape(string &filename, int *vshape, int *gshape, int *kshape) {
    ifstream fin(filename.c_str());
    string line;

    getline(fin, line, '\n');
    load_line(line, vshape, 2);

    getline(fin, line, '\n');
    load_line(line, gshape, 4);

    getline(fin, line, '\n');
    load_line(line, kshape, 5);
}


template<typename T>
void load_data(string filepath, T *arr, const int n) {
    ifstream fin(filepath.c_str(), ios::binary | ios::in);
    fin.read((char*)arr, sizeof(T) * n);
    fin.close();
}


template <typename T>
void print_arr(T *arr, int len) {
    if (sizeof(T) == 4 || sizeof(T) == 8) {
        for (int i = 0; i < len; i++) {
            cout << arr[i] << " ";
        }
    }
    else if (sizeof(T) == 16) {
        for (int i = 0; i < len; i++) {
            cout << creal(arr[i]) << (cimag(arr[i]) >= 0 ? "+" : "") << cimag(arr[i]) << "j ";
        }
    }
    cout << endl;
}


template<typename T>
double diff(T *a, T *b, int n) {
    double err = 0.0;
    for (int i = 0; i < n; i++) {
        err = max(err, cabs(a[i] - b[i]));
    }
    return err;
}


#endif