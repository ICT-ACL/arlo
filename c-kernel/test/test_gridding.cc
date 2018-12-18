// test/test_gridding.cc
// tests of convolutional gridding
// Author: You Haihang, Yang Runkai, Liu Tao,
// Applied Computing Lab, ICT, CAS

#include <iostream>
#include <string>
#include <ctime>

#include "convolutional_gridding.h"
#include "utils.h"

using namespace std;


void test_degrid(string &data_dir) {
    string shape_filepath = data_dir + "/shapes.txt";
    int vshape[2], gshape[4], kshape[5];
    load_shape(shape_filepath, vshape, gshape, kshape);
    int vlen = vshape[0] * vshape[1];
    int glen = gshape[0] * gshape[1] * gshape[2] * gshape[3];
    int klen = kshape[0] * kshape[1] * kshape[2] * kshape[3] * kshape[4];
    int nvis = vshape[0];

    complex_t *kernels = new complex_t[klen];
    int *kernel_indices = new int[nvis];
    complex_t *uvgrid = new complex_t[glen];
    double *vuvwmap = new double[nvis * 3];
    int *vfrequencymap = new int[nvis];
    complex_t *vis = new complex_t[vlen];
    complex_t *vis_exact = new complex_t[vlen];

    load_data(data_dir + "/" + "kernels.dat", kernels, klen);
    load_data(data_dir + "/" + "kernel_indices.dat", kernel_indices, nvis);
    load_data(data_dir + "/" + "uvgrid.dat", uvgrid, glen);
    load_data(data_dir + "/" + "vuvwmap.dat", vuvwmap, nvis * 3);
    load_data(data_dir + "/" + "vfrequencymap.dat", vfrequencymap, nvis);

    clock_t start = clock();
    convolutional_degrid(vis, vshape, kernels, kernel_indices, kshape, uvgrid, gshape,
                         vuvwmap, vfrequencymap);
    clock_t stop = clock();
    double wtime = 1.0 * (stop - start) / CLOCKS_PER_SEC;
    printf("Optimized Time: %.2lfs\n", wtime);

    load_data(data_dir + "/" + "vis.dat", vis_exact, vlen);
    double err = diff(vis, vis_exact, vlen);
    cout << "Degridding Error: " << err << endl;

    delete [] kernels;
    delete [] kernel_indices;
    delete [] uvgrid;
    delete [] vuvwmap;
    delete [] vfrequencymap;
    delete [] vis;
    delete [] vis_exact;
}


void test_grid(string &data_dir) {
    string shape_filepath = data_dir + "/shapes.txt";
    int vshape[2], gshape[4], kshape[5];
    load_shape(shape_filepath, vshape, gshape, kshape);
    int nvis = vshape[0], inchan = gshape[0], inpol = gshape[1];
    int vlen = vshape[0] * vshape[1];
    int glen = gshape[0] * gshape[1] * gshape[2] * gshape[3];
    int klen = kshape[0] * kshape[1] * kshape[2] * kshape[3] * kshape[4];
    int slen = inchan * inpol;

    complex_t *kernels = new complex_t[klen];
    int *kernel_indices = new int[nvis];
    complex_t *uvgrid = new complex_t[glen];
    complex_t *uvgrid_exact = new complex_t[glen];
    complex_t *vis = new complex_t[vlen];
    double *visweights = new double[vlen];
    double *vuvwmap = new double[nvis * 3];
    int *vfrequencymap = new int[nvis];
    double *sumwt = new double[slen];
    double *sumwt_exact = new double[slen];

    load_data(data_dir + "/" + "kernels.dat", kernels, klen);
    load_data(data_dir + "/" + "kernel_indices.dat", kernel_indices, nvis);
    load_data(data_dir + "/" + "uvgrid_before.dat", uvgrid, glen);
    load_data(data_dir + "/" + "vis.dat", vis, vlen);
    load_data(data_dir + "/" + "visweights.dat", visweights, vlen);
    load_data(data_dir + "/" + "vuvwmap.dat", vuvwmap, nvis * 3);
    load_data(data_dir + "/" + "vfrequencymap.dat", vfrequencymap, nvis);

    clock_t start = clock();
    convolutional_grid(uvgrid, gshape, sumwt, vis, vshape, visweights, kernels, kernel_indices, kshape,
                         vuvwmap, vfrequencymap);
    clock_t stop = clock();
    double wtime = 1.0 * (stop - start) / CLOCKS_PER_SEC;
    printf("Optimized Time: %.2lfs\n", wtime);

    load_data(data_dir + "/" + "uvgrid_after.dat", uvgrid_exact, glen);
    load_data(data_dir + "/" + "sumwt.dat", sumwt_exact, slen);

    double err_uvgrid = diff(uvgrid, uvgrid_exact, glen);
    double err_sumwt = diff(sumwt, sumwt_exact, slen);
    cout << "Gridding Error (uvgrid): " << err_uvgrid << endl;
    cout << "Gridding Error (sumwt) : " << err_sumwt << endl;

    delete [] kernels;
    delete [] kernel_indices;
    delete [] uvgrid;
    delete [] uvgrid_exact;
    delete [] vis;
    delete [] visweights;
    delete [] vuvwmap;
    delete [] vfrequencymap;
    delete [] sumwt;
    delete [] sumwt_exact;
}

int main(int argc, char *argv[]) {
    string data_dir = "./data";

    if (!strcmp(argv[1], "degrid")) {
        test_degrid(data_dir);
    }
    else if (!strcmp(argv[1], "grid")) {
        test_grid(data_dir);
    }

    return 0;
}