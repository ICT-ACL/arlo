// test/test_gridding.cc
// tests of convolutional gridding
// Author: You Haihang, Yang Runkai, Liu Tao,
// Applied Computing Lab, ICT, CAS

#include <iostream>
#include <string>
#include <ctime>
#include <getopt.h>

#include "convolutional_gridding.h"
#include "utils.h"

using namespace std;


void test_degrid(string &data_dir, vector<int> &vshape, vector<int> &gshape, vector<int> &kshape) {
    int vlen = prod(vshape), glen = prod(gshape), klen = prod(kshape);
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
    convolutional_degrid(vis, &vshape[0], kernels, kernel_indices, &kshape[0], uvgrid, &gshape[0],
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


void test_grid(string &data_dir, vector<int> &vshape, vector<int> &gshape, vector<int> &kshape) {
    int vlen = prod(vshape), glen = prod(gshape), klen = prod(kshape);
    int nvis = vshape[0], inchan = gshape[0], inpol = gshape[1];
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
    convolutional_grid(uvgrid, &gshape[0], sumwt, vis, &vshape[0], visweights, kernels, kernel_indices, &kshape[0],
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
    string mode = "grid";
    string vshape_str = "0,0";
    string gshape_str = "0,0,0,0";
    string kshape_str = "0,0,0,0,0";

    int c;
    while (true) {
        int option_index = 0;
        static struct option long_options[] = {
            {"data_dir", required_argument, 0, 1},
            {"mode", required_argument, 0, 2},
            {"vshape", required_argument, 0, 3},
            {"gshape", required_argument, 0, 4},
            {"kshape", required_argument, 0, 5}
        };

        c = getopt_long(argc, argv, "", long_options, &option_index);
        if (c == -1) break;

        switch (c) {
            case 1: data_dir = string(optarg); break;
            case 2: mode = string(optarg); break;
            case 3: vshape_str = string(optarg); break;
            case 4: gshape_str = string(optarg); break;
            case 5: kshape_str = string(optarg); break;
            default: cout << "Wrong options." << endl; break;
        }
    }

    vector<int> vshape = str2arr(vshape_str);
    vector<int> gshape = str2arr(gshape_str);
    vector<int> kshape = str2arr(kshape_str);

    if (vshape.size() == 2 && gshape.size() == 4 && kshape.size() == 5) {
        if (mode == "degrid") {
            test_degrid(data_dir, vshape, gshape, kshape);
        }
        else if (mode == "grid") {
            test_grid(data_dir, vshape, gshape, kshape);
        }
        else {}
    }

    return 0;
}