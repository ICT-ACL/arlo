#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
// #include <cblas.h>

#include "convolutional_gridding.h"
#include "debug.h"


using namespace std;

/**
* Compute whole and fractional parts of coordinates, rounded to
* kernel_oversampling-th fraction of pixel size
* The fractional values are rounded to nearest 1/kernel_oversampling pixel value. At
* fractional values greater than (kernel_oversampling-0.5)/kernel_oversampling coordinates are
* rounded to next integer index.
* :param npixel: Number of pixels in total
* :param kernel_oversampling: Fractional values to round to
* :param p: Coordinate in range [-.5,.5)
*/
void frac_coord(int *flx, int *fracx, double *p, const int length, const int npixel, int kernel_oversampling) {
	for (int i = 0; i < length; i++) {
		assert(p[i] >= -0.5 && p[i] < 0.5);
	}

	double *x = (double*)malloc(sizeof(double) * length);
	for (int i = 0; i < length; i++) {
		x[i] = npixel / 2 + npixel * p[i];
		flx[i] = floor(x[i] + 0.5 / kernel_oversampling);
		fracx[i] = round((x[i] - flx[i]) * kernel_oversampling);
	}

	free(x);
}


// calculate coordinates
void cal_coord(int *a, int *af, double *vuvwmap, 
			   const int length, const int npixel, const int axis, 
			   int kernel_oversampling, const int g) {

	double *map = (double*)malloc(sizeof(double) * length);

	for (int i = 0; i < length; i++) {
		map[i] = vuvwmap[i * 3 + axis];
	}

	frac_coord(a, af, map, length, npixel, kernel_oversampling);


	for (int i = 0; i < length; i++) {
		a[i] -= g / 2;
	}

	free(map);
}



/**
* Convolutional degridding with frequency and polarisation independent
* Takes into account fractional `uv` coordinate values where the GCF
* is oversampled
* :param kernel_indices: indices of convolution kernel
* :param kernels: list of oversampled convolution kernel
* :param uvgrid:   The uv plane to de-grid from
* :param kshape: Shape of kernels
* :param vshape: Shape of visibility
* :param gshape: Shape of uvgrid
* :param vuvwmap: function to map uvw to grid fractions
* :param vfrequencymap: function to map frequency to image channels
* :param vpolarisationmap: function to map polarisation to image polarisation
* :return: Array of visibilities.
*/
void convolutional_degrid(complex_t *vis, const int *vshape, 
					      complex_t* kernels, int *kernel_indices, const int *kshape, 
						  complex_t *uvgrid, const int *gshape, 
						  double *vuvwmap, int *vfrequencymap) {

	int nvis = vshape[0], vnpol = vshape[1];
	int nkernels = kshape[0], kernel_oversampling = kshape[1], gh = kshape[3], gw = kshape[4];
	int ksize = kernel_oversampling * kernel_oversampling * gh * gw;
	int ny = gshape[2], nx = gshape[3];

	assert(gh % 2 == 0);
	assert(gw % 2 == 0);

	// uvw -> fraction of grid mapping
	int *y  = (int*)malloc(sizeof(int) * nvis);
	int *yf = (int*)malloc(sizeof(int) * nvis);
	int *x  = (int*)malloc(sizeof(int) * nvis);
	int *xf = (int*)malloc(sizeof(int) * nvis);
	cal_coord(y, yf, vuvwmap, nvis, ny, 1, kernel_oversampling, gh);
	cal_coord(x, xf, vuvwmap, nvis, nx, 0, kernel_oversampling, gw);

	// conjugate of kernels
	complex_t *ckernels = (complex_t*)malloc(sizeof(complex_t) * ksize * nkernels);
	for (int i = 0; i < ksize * nkernels; i++) {
		ckernels[i] = conj(kernels[i]);
	}

	memset(vis, 0, sizeof(complex_t) * vnpol * nvis);
	
	// degridding
	// for (int pol = 0; pol < vnpol; pol++) {

	// 	for (int s = 0; s < nvis; s++) {

	// 		int kind = kernel_indices[s];
	// 		int chan = vfrequencymap[s];
	// 		int yy = y[s], yyf = yf[s];
	// 		int xx = x[s], xxf = xf[s];

	// 		complex_t* grid = uvgrid + (chan * vnpol + pol) * ny * nx;
	// 		complex_t* kernel = ckernels + kind * ksize + (yyf * kernel_oversampling + xxf) * gh * gw;

	// 		for (int i = 0; i < gh; i++) {
	// 			for (int j = 0; j < gw; j++) {

	// 				int v = yy + i, u = xx + j;
	// 				vis[s*vnpol+pol] += grid[v * nx + u] * kernel[i * gw + j];
	// 			}
	// 		}
	// 	}
	// }

	#pragma omp parallel for
	for (int ps = 0; ps < vnpol * nvis; ps++) {
		int pol = ps / nvis;
		int s = ps % nvis;

		int kind = kernel_indices[s];
		int chan = vfrequencymap[s];
		int yy = y[s], yyf = yf[s];
		int xx = x[s], xxf = xf[s];

		complex_t* grid = uvgrid + (chan * vnpol + pol) * ny * nx;
		complex_t* kernel = ckernels + kind * ksize + (yyf * kernel_oversampling + xxf) * gh * gw;

		for (int i = 0; i < gh; i++) {
			for (int j = 0; j < gw; j++) {

				int v = yy + i, u = xx + j;
				vis[s*vnpol+pol] += grid[v * nx + u] * kernel[i * gw + j];
			}
		}
	}

	free(y);
	free(yf);
	free(x);
	free(xf);
	free(ckernels);
}



/**
* Grid after convolving with frequency and polarisation independent gcf
* Takes into account fractional `uv` coordinate values where the GCF is oversampled
* :param kernels: List of oversampled convolution kernels
* :param uvgrid: Grid to add to [nchan, npol, npixel, npixel]
* :param vis: Visibility values
* :param visweights: Visibility weights
* :param vuvwmap: map uvw to grid fractions
* :param vfrequencymap: map frequency to image channels
* :param vpolarisationmap: map polarisation to image polarisation
* :return: uv grid[nchan, npol, ny, nx], sumwt[nchan, npol]
*/
void convolutional_grid(complex_t *uvgrid, const int *gshape, 
					    double *sumwt,
					    complex_t *vis, const int *vshape, 
					    double *visweights,
				        complex_t* kernels, int *kernel_indices, const int *kshape, 
					    double *vuvwmap, int *vfrequencymap) {
	
	int nvis = vshape[0], vnpol = vshape[1];
	int nkernels = kshape[0], kernel_oversampling = kshape[1], gh = kshape[3], gw = kshape[4];
	int ksize = kernel_oversampling * kernel_oversampling * gh * gw;
	int inchan = gshape[0], inpol = gshape[1], ny = gshape[2], nx = gshape[3];

	assert(gh % 2 == 0);
	assert(gw % 2 == 0);

	// uvw -> fraction of grid mapping
	int *y  = (int*)malloc(sizeof(int) * nvis);
	int *yf = (int*)malloc(sizeof(int) * nvis);
	int *x  = (int*)malloc(sizeof(int) * nvis);
	int *xf = (int*)malloc(sizeof(int) * nvis);
	cal_coord(y, yf, vuvwmap, nvis, ny, 1, kernel_oversampling, gh);
	cal_coord(x, xf, vuvwmap, nvis, nx, 0, kernel_oversampling, gw);

	double *wts = visweights;
	complex_t *viswt = (complex_t*)malloc(sizeof(complex_t) * nvis * vnpol);
	for (int i = 0; i < vnpol * nvis; i++) {
		viswt[i] = vis[i] * visweights[i];
	}

	memset(uvgrid, 0, sizeof(complex_t) * inchan * inpol * ny * nx);
	memset(sumwt, 0, sizeof(double) * inchan * inpol);

	// gridding
	// for (int pol = 0; pol < vnpol; pol++) {

	// 	for (int s = 0; s < nvis; s++) {

	// 		int kind = kernel_indices[s];
	// 		int chan = vfrequencymap[s];
	// 		int yy = y[s], yyf = yf[s];
	// 		int xx = x[s], xxf = xf[s];

	// 		sumwt[chan * vnpol + pol] += wts[s * vnpol + pol];

	// 		complex_t* grid = uvgrid + (chan * vnpol + pol) * ny * nx;
	// 		complex_t* kernel = kernels + kind * ksize + (yyf * kernel_oversampling + xxf) * gh * gw;
	// 		complex_t vt = viswt[s * vnpol + pol];

	// 		for (int i = 0; i < gh; i++) {
	// 			for (int j = 0; j < gw; j++) {

	// 				int v = yy + i, u = xx + j;
	// 				grid[v * nx + u] += kernel[i * gw + j] * vt;

	// 			}
	// 		}
	// 	}
	// }


	for (int ps = 0; ps < vnpol * nvis; ps++) {
		int pol = ps / nvis;
		int s = ps % nvis;

		int kind = kernel_indices[s];
		int chan = vfrequencymap[s];
		int yy = y[s], yyf = yf[s];
		int xx = x[s], xxf = xf[s];

		sumwt[chan * vnpol + pol] += wts[s * vnpol + pol];

		complex_t* grid = uvgrid + (chan * vnpol + pol) * ny * nx;
		complex_t* kernel = kernels + kind * ksize + (yyf * kernel_oversampling + xxf) * gh * gw;
		complex_t vt = viswt[s * vnpol + pol];


		for (int i = 0; i < gh; i++) {
			// #pragma simd
			for (int j = 0; j < gw; j++) {

				int v = yy + i, u = xx + j;
				grid[v * nx + u] += kernel[i * gw + j] * vt;

			}
		}

		// for (int i = 0; i < gh; i++) {
		// 	cblas_zaxpy(gw, &vt, kernel + i*gw, 1, grid + ((yy+i)*nx + xx), 1);
		// }
	}
	
	free(y);
	free(yf);
	free(x);
	free(xf);
}


int convolutional_degrid_c(PyArrayObject *&vis_obj,
                 PyArrayObject *&kernels_obj,
                 PyArrayObject *&kernel_indices_obj,
                 PyArrayObject *&uvgrid_obj,
				 PyArrayObject *&vuvwmap_obj,
				 PyArrayObject *&vfrequencymap_obj) {

    complex_t *vis = (complex_t*)vis_obj->data;
	complex_t *kernels = (complex_t*)kernels_obj->data;
	int *kernel_indices = (int*)kernel_indices_obj->data;
	complex_t *uvgrid = (complex_t*)uvgrid_obj->data;
	double *vuvwmap = (double*)vuvwmap_obj->data;
	int *vfrequencymap = (int*)vfrequencymap_obj->data;

	int *vshape = (int*)malloc(sizeof(int) * vis_obj->nd);
	int *kshape = (int*)malloc(sizeof(int) * kernels_obj->nd);
	int *gshape = (int*)malloc(sizeof(int) * uvgrid_obj->nd);

	for (int i = 0; i < vis_obj->nd; i++) {
		vshape[i] = vis_obj->dimensions[i];
	}
	for (int i = 0; i < kernels_obj->nd; i++) {
		kshape[i] = kernels_obj->dimensions[i];
	}
	for (int i = 0; i < uvgrid_obj->nd; i++) {
		gshape[i] = uvgrid_obj->dimensions[i];
	}

	convolutional_degrid(vis, vshape, 
						kernels, kernel_indices, kshape, 
						uvgrid, gshape, 
						vuvwmap, vfrequencymap);

	return 0;
}


int convolutional_grid_c(PyArrayObject *&uvgrid_obj,
						 PyArrayObject *&sumwt_obj,
						 PyArrayObject *&vis_obj,
						 PyArrayObject *&visweights_obj,
		                 PyArrayObject *&kernels_obj,
		                 PyArrayObject *&kernel_indices_obj,
						 PyArrayObject *&vuvwmap_obj,
						 PyArrayObject *&vfrequencymap_obj) {

	complex_t *uvgrid = (complex_t*)uvgrid_obj->data;
	double *sumwt = (double*)sumwt_obj->data;
	complex_t *vis = (complex_t*)vis_obj->data;
	double *visweights = (double*)visweights_obj->data;
	complex_t *kernels = (complex_t*)kernels_obj->data;
	int *kernel_indices = (int*)kernel_indices_obj->data;
	double *vuvwmap = (double*)vuvwmap_obj->data;
	int *vfrequencymap = (int*)vfrequencymap_obj->data;

	int *vshape = (int*)malloc(sizeof(int) * vis_obj->nd);
	int *kshape = (int*)malloc(sizeof(int) * kernels_obj->nd);
	int *gshape = (int*)malloc(sizeof(int) * uvgrid_obj->nd);

	for (int i = 0; i < vis_obj->nd; i++) {
		vshape[i] = vis_obj->dimensions[i];
	}
	for (int i = 0; i < kernels_obj->nd; i++) {
		kshape[i] = kernels_obj->dimensions[i];
	}
	for (int i = 0; i < uvgrid_obj->nd; i++) {
		gshape[i] = uvgrid_obj->dimensions[i];
	}

	convolutional_grid(uvgrid, gshape, sumwt,
					    vis, vshape, visweights,
				        kernels, kernel_indices, kshape, 
					    vuvwmap, vfrequencymap);

	return 0;
}