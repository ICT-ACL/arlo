#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <string>
#include <omp.h>

#include "cleaners.h"

#define MAX_LEN 1024

using namespace std;

template<typename T>
void load_data(const char *data_dir, const char *data_name, T *arr, const int len) {
	char filename[MAX_LEN];
	strcat(strcpy(filename, data_dir), data_name);

	FILE *fr = fopen(filename, "r");
	if (fr == NULL) {
		fprintf(stderr, "no such file: %s\n", filename);
		exit(-1);
	}

	if (sizeof(T) == sizeof(int)) {

		long *tmp = (long*)malloc(sizeof(long) * len);
		size_t status = fread(tmp, sizeof(long), len, fr);
		for (int i = 0; i < len; i++) {
			arr[i] = int(tmp[i]);
		}
		free(tmp);
	}
	else {
		size_t status = fread(arr, sizeof(T), len, fr);
	}

	fclose(fr);
}

template<typename T>
void store_data(const char *data_dir, const char *data_name, T *arr, const int len) {
	char filename[MAX_LEN];
	strcat(strcpy(filename, data_dir), data_name);

	FILE *fw = fopen(filename, "w");

	if (sizeof(T) == sizeof(int)) {
		long *tmp = (long*)malloc(sizeof(long) * len);
		for (int i = 0; i < len; i++) {
			tmp[i] = long(arr[i]);
		}
		fwrite(tmp, sizeof(long), len, fw);
		free(tmp);
	}
	else {
		fwrite(arr, sizeof(T), len, fw);
	}
	fclose(fw);
}

double check(double *a, double *b, int length) {
	double max = 0.0;

	for (int i = 0; i < length; i++) {
		double val = fabs(a[i] - b[i]);
		max = max > val ? max : val;
	}

	return max;
}


int main(int argc, char *argv[]) {
	const char data_root_dir[MAX_LEN] = "../data/";
	char data_dir[MAX_LEN];
	strcat(strcpy(data_dir, data_root_dir), "clean/");

	double gain = 0.7;
	double absolutethresh = 0.453597533897;
	string findpeak = "ARL";
	int niter = 1;
	if (argc = 2) {
		niter = atoi(argv[1]);
	}

	const int nscales = 4;
	const int nmoments = 3;
	const int nx = 512;
	const int ny = 512;

	// output
	double *m_model = (double*)malloc(sizeof(double) * nmoments * ny * nx);
	double *m_model_exact = (double*)malloc(sizeof(double) * nmoments * ny * nx);
	double *residual = (double*)malloc(sizeof(double) * nmoments * ny * nx);
	double *residual_exact = (double*)malloc(sizeof(double) * nmoments * ny * nx);
	// input
	double *ldirty = (double*)malloc(sizeof(double) * nmoments * ny * nx);
	double *psf = (double*)malloc(sizeof(double) * 2*nmoments * ny * nx);
	double *scalestack = (double*)malloc(sizeof(double) * nscales * ny * nx);
	double *ssmmpsf = (double*)malloc(sizeof(double) * nscales * nscales * nmoments * nmoments * ny * nx);
	double *hsmmpsf = (double*)malloc(sizeof(double) * nscales * nmoments * nmoments);
	double *ihsmmpsf = (double*)malloc(sizeof(double) * nscales * nmoments * nmoments);
	double *smresidual = (double*)malloc(sizeof(double) * nscales * nmoments * ny * nx);
	double *windowstack = NULL;

	load_data(data_dir, "m_model.dat", m_model_exact, nmoments * ny * nx);
	load_data(data_dir, "residual.dat", residual_exact, nmoments * ny * nx);
	load_data(data_dir, "ldirty.dat", ldirty, nmoments * ny * nx);
	load_data(data_dir, "psf.dat", psf, 2*nmoments * ny * nx);
	load_data(data_dir, "scalestack.dat", scalestack, nscales * ny * nx);
	load_data(data_dir, "ssmmpsf.dat", ssmmpsf, nscales * nscales * nmoments * nmoments * ny * nx);
	load_data(data_dir, "hsmmpsf.dat", hsmmpsf, nscales * nmoments * nmoments);
	load_data(data_dir, "ihsmmpsf.dat", ihsmmpsf, nscales * nmoments * nmoments);
	load_data(data_dir, "smresidual.dat", smresidual, nscales * nmoments * ny * nx);

	double time = -omp_get_wtime();
	msmfsclean_kernel(m_model, residual, scalestack, smresidual, ssmmpsf, hsmmpsf, ihsmmpsf, ldirty, psf, 
					  nscales, nmoments, nx, ny, windowstack, gain, absolutethresh, niter, findpeak);
	time += omp_get_wtime();

	double err_clean = check(m_model, m_model_exact, nmoments * ny * nx);
	double err_residual = check(residual, residual_exact, nmoments * ny * nx);
	printf("error_clean: %lg\n", err_clean);
	printf("error_residual: %lg\n", err_residual);
	printf("time: %.4lfs\n", time);

	return 0;
}