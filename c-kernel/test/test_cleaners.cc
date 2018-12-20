// test/test_clean.cc
// tests of cleaners
// Author: You Haihang, Yang Runkai, Liu Tao,
// Applied Computing Lab, ICT, CAS 

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <getopt.h>
#include <string>

#include "cleaners.h"
#include "utils.h"

using namespace std;


void test_cleaners(string &data_dir, int niter, double gain, double thresh,
                   double fracthresh, const char *findpeak,
                   const int nscales, const int nmoments,
                   const int nx, const int ny) {

	double *m_model = new double[nmoments * ny * nx];
	double *m_model_exact = new double[nmoments * ny * nx];
	double *residual = new double[nmoments * ny * nx];
	double *residual_exact = new double[nmoments * ny * nx];
	double *ldirty = new double[nmoments * ny * nx];
	double *psf = new double[2*nmoments * ny * nx];
	double *pscalestack = new double[nscales * ny * nx];
	double *ssmmpsf = new double[nscales * nscales * nmoments * nmoments * ny * nx];
	double *hsmmpsf = new double[nscales * nmoments * nmoments];
	double *ihsmmpsf = new double[nscales * nmoments * nmoments];
	double *smresidual = new double[nscales * nmoments * ny * nx];
	double *windowstack = NULL;

	load_data(data_dir + "/" + "m_model.dat", m_model_exact, nmoments * ny * nx);
	load_data(data_dir + "/" + "residual.dat", residual_exact, nmoments * ny * nx);
	load_data(data_dir + "/" + "ldirty.dat", ldirty, nmoments * ny * nx);
	load_data(data_dir + "/" + "psf.dat", psf, 2*nmoments * ny * nx);
	load_data(data_dir + "/" + "pscalestack.dat", pscalestack, nscales * ny * nx);
	load_data(data_dir + "/" + "ssmmpsf.dat", ssmmpsf, nscales * nscales * nmoments * nmoments * ny * nx);
	load_data(data_dir + "/" + "hsmmpsf.dat", hsmmpsf, nscales * nmoments * nmoments);
	load_data(data_dir + "/" + "ihsmmpsf.dat", ihsmmpsf, nscales * nmoments * nmoments);
	load_data(data_dir + "/" + "smresidual.dat", smresidual, nscales * nmoments * ny * nx);

	double absolutethresh = 0.0;
    for (int i = 0; i < nx * ny; i++) {
        absolutethresh = max(absolutethresh, fabs(smresidual[i]));
    }
    absolutethresh = max(thresh, fracthresh * absolutethresh);

	clock_t start = clock();
	msmfsclean_kernel(m_model, residual, pscalestack, smresidual, ssmmpsf, hsmmpsf, ihsmmpsf, ldirty, psf,
					  nscales, nmoments, nx, ny, windowstack, gain, absolutethresh, niter, findpeak);
	clock_t stop = clock();
	double wtime = 1.0 * (stop - start) / CLOCKS_PER_SEC;
	printf("Optimized Time: %.2lfs\n", wtime);

	double err_model = diff_relative(m_model, m_model_exact, nmoments * ny * nx);
	double err_residual = diff(residual, residual_exact, nmoments * ny * nx);
	cout << "Cleaners Error (m_model, relative) : " << err_model << endl;
	cout << "Cleaners Error (residual): " << err_residual << endl;

	delete [] m_model;
    delete [] m_model_exact;
    delete [] residual;
    delete [] residual_exact;
    delete [] ldirty;
    delete [] psf;
    delete [] pscalestack;
    delete [] ssmmpsf;
    delete [] hsmmpsf;
    delete [] ihsmmpsf;
    delete [] smresidual;
}


int main(int argc, char *argv[]) {
    string data_dir = "./data";
    int niter = 0;
	double gain = 0.0;
	double thresh=0.0;
	double fracthresh=0.0;
	int nscales = 0;
	int nmoments = 0;
	int nx = 0;
	int ny = 0;

    int c;
    while (true) {
        int option_index = 0;
        static struct option long_options[] = {
            {"data_dir", required_argument, 0, 1},
            {"niter", required_argument, 0, 2},
            {"gain", required_argument, 0, 3},
            {"thresh", required_argument, 0, 4},
            {"fracthresh", required_argument, 0, 5},
            {"nscales", required_argument, 0, 6},
            {"nmoments", required_argument, 0, 7},
            {"nx", required_argument, 0, 8},
            {"ny", required_argument, 0, 9}
        };

        c = getopt_long(argc, argv, "", long_options, &option_index);

        if (c == -1) break;

        switch (c) {
            case 1: data_dir = string(optarg); break;
            case 2: niter = atoi(optarg); break;
            case 3: gain = atof(optarg); break;
            case 4: thresh = atof(optarg); break;
            case 5: fracthresh = atof(optarg); break;
            case 6: nscales = atof(optarg); break;
            case 7: nmoments = atof(optarg); break;
            case 8: nx = atof(optarg); break;
            case 9: ny = atof(optarg); break;
            default: cout << "Wrong options" << endl; break;
        }
    }

    test_cleaners(data_dir, niter, gain, thresh, fracthresh, "ARL", nscales, nmoments, nx, ny);

    return 0;
}