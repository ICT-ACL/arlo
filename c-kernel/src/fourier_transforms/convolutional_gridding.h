// fourier_transforms/convolutional_gridding.h
// convolutional gridding & degridding codes
// Author: You Haihang, Yang Runkai, Liu Tao,
// Applied Computing Lab, ICT, CAS 
#ifndef CONVOLUTIONAL_GRIDDING_H
#define CONVOLUTIONAL_GRIDDING_H

#include <complex.h>
#include <ndarrayobject.h>

typedef double complex complex_t;

// #include "../utils/fftw.h"

void convolutional_degrid(complex_t *vis, const int *vshape, 
					      complex_t* kernels, int *kernel_indices, const int *kshape, 
						  complex_t *uvgrid, const int *gshape, 
						  double *vuvwmap, int *vfrequencymap);

void convolutional_grid(complex_t *uvgrid, const int *gshape, 
					    double *sumwt,
					    complex_t *vis, const int *vshape, 
					    double *visweights,
				        complex_t* kernels, int *kernel_indices, const int *kshape, 
					    double *vuvwmap, int *vfrequencymap);

int convolutional_degrid_c(PyArrayObject *&vis_obj,
                 PyArrayObject *&kernels_obj,
                 PyArrayObject *&kernel_indices_obj,
                 PyArrayObject *&uvgrid_obj,
				 PyArrayObject *&vuvwmap_obj,
				 PyArrayObject *&vfrequencymap_obj);

int convolutional_grid_c(PyArrayObject *&uvgrid_obj,
						 PyArrayObject *&sumwt_obj,
						 PyArrayObject *&vis_obj,
						 PyArrayObject *&visweights_obj,
		                 PyArrayObject *&kernels_obj,
		                 PyArrayObject *&kernel_indices_obj,
						 PyArrayObject *&vuvwmap_obj,
						 PyArrayObject *&vfrequencymap_obj);				 
#endif
