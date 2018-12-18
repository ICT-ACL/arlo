// test/test_clean.cc
// tests of cleaners
// Author: You Haihang, Yang Runkai, Liu Tao,
// Applied Computing Lab, ICT, CAS 

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <string>
#include <omp.h>

#include "cleaners.h"
#include "utils.h"

using namespace std;


void test_cleaners(string &data_dir) {
	double gain = 0.7;
	double absolutethresh = 0.19999971650614415;
	const char findpeak[] = "ARL";
	int niter = 2;

	const int nscales = 4;
	const int nmoments = 3;
	const int nx = 512;
	const int ny = 512;

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

	clock_t start = clock();
	msmfsclean_kernel(m_model, residual, pscalestack, smresidual, ssmmpsf, hsmmpsf, ihsmmpsf, ldirty, psf,
					  nscales, nmoments, nx, ny, windowstack, gain, absolutethresh, niter, findpeak);
	clock_t stop = clock();
	double wtime = 1.0 * (stop - start) / CLOCKS_PER_SEC;
	printf("Optimized Time: %.2lfs\n", wtime);

	double err_model = diff(m_model, m_model_exact, nmoments * ny * nx);
	double err_residual = diff(residual, residual_exact, nmoments * ny * nx);
	cout << "Cleaners Error (m_model) : " << err_model << endl;
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
    test_cleaners(data_dir);

    return 0;
}