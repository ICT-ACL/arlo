#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>

#include "convolutional_gridding.h"

#define MAX_LEN 1024


void load_shape(const char *data_dir, const char *data_name, 
				int *vshape, int *gshape, int *kshape,
				const int n_vshapes, const int n_gshapes, const int n_kshapes) {

	char filename[MAX_LEN];
	strcat(strcpy(filename, data_dir), data_name);

	FILE *fr = fopen(filename, "r");
	if (fr == NULL) {
		fprintf(stderr, "no such file: %s\n", filename);
		exit(-1);
	}
	for (int i = 0; i < n_vshapes; i++) {
		fscanf(fr, "%d", &vshape[i]);
	}
	for (int i = 0; i < n_gshapes; i++) {
		fscanf(fr, "%d", &gshape[i]);
	}
	for (int i = 0; i < n_kshapes; i++) {
		fscanf(fr, "%d", &kshape[i]);
	}

	fclose(fr);
}


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
		fread(tmp, sizeof(long), len, fr);
		for (int i = 0; i < len; i++) {
			arr[i] = int(tmp[i]);
		}
		free(tmp);
	}
	else {
		fread(arr, sizeof(T), len, fr);
	}

	fclose(fr);
}


void load_data(const char *data_dir, const char *data_name, complex_t *arr, const int len) {
	char filename[MAX_LEN];
	strcat(strcpy(filename, data_dir), data_name);

	FILE *fr = fopen(filename, "r");
	if (fr == NULL) {
		fprintf(stderr, "no such file: %s\n", filename);
		exit(-1);
	}

	fread(arr, sizeof(complex_t), len, fr);

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

void store_data(const char *data_dir, const char *data_name, complex_t *arr, const int len) {
	char filename[MAX_LEN];
	strcat(strcpy(filename, data_dir), data_name);

	FILE *fw = fopen(filename, "w");
	fwrite(arr, sizeof(complex_t), len, fw);
	fclose(fw);
}



double check(complex_t *a, complex_t *b, int length) {
	double max = 0.0;

	for (int i = 0; i < length; i++) {
		double val = cabs(a[i] - b[i]);
		max = max > val ? max : val;
	}

	return max;
}

double check(double *a, double *b, int length) {
	double max = 0.0;

	for (int i = 0; i < length; i++) {
		double val = fabs(a[i] - b[i]);
		max = max > val ? max : val;
	}

	return max;
}


double test_convolutional_degrid(const char *data_root_dir, double &err_degrid) {
	// data directory
	char data_dir[MAX_LEN];
	strcat(strcpy(data_dir, data_root_dir), "degrid/");

	// load data shapes
	const int n_vshapes = 2, n_gshapes = 4, n_kshapes = 5;
	int vshape[n_vshapes], gshape[n_gshapes], kshape[n_kshapes];
	load_shape(data_dir, "shape.dat", vshape, gshape, kshape, n_vshapes, n_gshapes, n_kshapes);

	int glength = gshape[0] * gshape[1] * gshape[2] * gshape[3];
	int klength = kshape[0] * kshape[1] * kshape[2] * kshape[3] * kshape[4];
	int vlength = vshape[0] * vshape[1];
	int nvis = vshape[0];

	complex_t *vis = (complex_t*)malloc(sizeof(complex_t) * vlength);
	complex_t *vis_exact = (complex_t*)malloc(sizeof(complex_t) * vlength);
	complex_t *uvgrid = (complex_t*)malloc(sizeof(complex_t) * glength);
	complex_t *kernels = (complex_t*)malloc(sizeof(complex_t) * klength);
	int *kernel_indices = (int*)malloc(sizeof(int) * nvis);
	int *vfrequencymap = (int*)malloc(sizeof(int) * nvis);
	double *vuvwmap = (double*)malloc(sizeof(double) * nvis * 3);

	// load data
	load_data(data_dir, "uvgrid.dat", uvgrid, glength);
	load_data(data_dir, "kernels.dat", kernels, klength);
	load_data(data_dir, "vis.dat", vis_exact, vlength);
	load_data(data_dir, "kernel_indices.dat", kernel_indices, nvis);
	load_data(data_dir, "vfrequencymap.dat", vfrequencymap, nvis);
	load_data(data_dir, "vuvwmap.dat", vuvwmap, nvis * 3);

	printf("grid:(%d,%d,%d,%d) = %d\n", gshape[0], gshape[1], gshape[2], gshape[3], glength);
	printf("kernels:(%d,%d,%d,%d,%d) = %d\n", kshape[0], kshape[1], kshape[2], kshape[3], kshape[4], klength);
	printf("vis:(%d,%d) = %d\n", vshape[0], vshape[1], vlength);

	// degridding
	double time = -omp_get_wtime();
	convolutional_degrid(vis, vshape, kernels, kernel_indices, kshape, 
						  uvgrid, gshape, vuvwmap, vfrequencymap);
	time += omp_get_wtime();

	// compute error (norm-inf)
	err_degrid = check(vis, vis_exact, vlength);

	free(vis);
	free(vis_exact);
	free(uvgrid);
	free(kernels);
	free(kernel_indices);
	free(vfrequencymap);
	free(vuvwmap);

	return time;
}


double test_convolutional_grid(const char *data_root_dir, double &err_grid, double &err_sum) {
	// data directory
	char data_dir[MAX_LEN];
	strcat(strcpy(data_dir, data_root_dir), "grid/");

	// load data shapes
	const int n_vshapes = 2, n_gshapes = 4, n_kshapes = 5;
	int vshape[n_vshapes], gshape[n_gshapes], kshape[n_kshapes];
	load_shape(data_dir, "shape.dat", vshape, gshape, kshape, n_vshapes, n_gshapes, n_kshapes);

	int glength = gshape[0] * gshape[1] * gshape[2] * gshape[3];
	int klength = kshape[0] * kshape[1] * kshape[2] * kshape[3] * kshape[4];
	int vlength = vshape[0] * vshape[1];
	int nfreq = gshape[0], nvis = vshape[0], nvpol = vshape[1];

	complex_t *uvgrid = (complex_t*)malloc(sizeof(complex_t) * glength);
	complex_t *uvgrid_exact = (complex_t*)malloc(sizeof(complex_t) * glength);
	double *sumwt = (double*)malloc(sizeof(double) * nfreq * nvpol);
	double *sumwt_exact = (double*)malloc(sizeof(double) * nfreq * nvpol);
	complex_t *kernels = (complex_t*)malloc(sizeof(complex_t) * klength);
	int *kernel_indices = (int*)malloc(sizeof(int) * nvis);
	complex_t *vis = (complex_t*)malloc(sizeof(complex_t) * vlength);
	double *vuvwmap = (double*)malloc(sizeof(double) * nvis * 3);
	int *vfrequencymap = (int*)malloc(sizeof(int) * nvis);
	double *visweights = (double*)malloc(sizeof(double) * nvis * nvpol);

	// load data
	load_data(data_dir, "uvgrid.dat", uvgrid_exact, glength);
	load_data(data_dir, "sumwt.dat", sumwt_exact, nfreq * nvpol);
	load_data(data_dir, "kernels.dat", kernels, klength);
	load_data(data_dir, "kernel_indices.dat", kernel_indices, nvis);
	load_data(data_dir, "vis.dat", vis, vlength);
	load_data(data_dir, "vuvwmap.dat", vuvwmap, nvis * 3);
	load_data(data_dir, "vfrequencymap.dat", vfrequencymap, nvis);
	load_data(data_dir, "visweights.dat", visweights, nvis * nvpol);

	printf("grid:(%d,%d,%d,%d) = %d\n", gshape[0], gshape[1], gshape[2], gshape[3], glength);
	printf("kernels:(%d,%d,%d,%d,%d) = %d\n", kshape[0], kshape[1], kshape[2], kshape[3], kshape[4], klength);
	printf("vis:(%d,%d) = %d\n", vshape[0], vshape[1], vlength);

	// gridding
	double time = -omp_get_wtime();
	convolutional_grid(uvgrid, gshape, sumwt, vis, vshape, visweights, 
						 kernels, kernel_indices, kshape, vuvwmap, vfrequencymap);
	time += omp_get_wtime();

	// compute error(norm-inf)
	err_grid = check(uvgrid, uvgrid_exact, glength);
	err_sum = check(sumwt, sumwt_exact, nfreq * nvpol);

	free(uvgrid);
	free(uvgrid_exact);
	free(sumwt);
	free(sumwt_exact);
	free(kernels);
	free(kernel_indices);
	free(vis);
	free(vuvwmap);
	free(vfrequencymap);
	free(visweights);

	return time;
}


int main(int argc, char *argv[]) {
	const char data_root_dir[MAX_LEN] = "../data/";

	int mode = 0;
	if (argc >= 2) {
		mode = atoi(argv[1]);
	}

	if (mode == 0) {	// degridding
		printf("degridding...\n");
		double err_degrid;
		double time_degrid = test_convolutional_degrid(data_root_dir, err_degrid);
		printf("degridding time: %.3lfs\n", time_degrid);
		printf("vis error:       %lg\n", err_degrid);
	}
	else {				// gridding
		printf("gridding...\n");
		double err_grid, err_sum;
		double time_grid = test_convolutional_grid(data_root_dir, err_grid, err_sum);
		printf("gridding time: %.3lfs\n", time_grid);
		printf("uvgrid error:  %lg\n", err_grid);
		printf("sumwt error:   %lg\n", err_sum);
	}
	
	return 0;
}