#ifndef CLEANERS_H
#define CLEANERS_H

#include <string>
#include <ndarrayobject.h>

using namespace std;

void msmfsclean_kernel(double *m_model, double *residual,
					   double *scalestack, double *smresidual,
					   double *ssmmpsf, double *hsmmpsf, double *ihsmmpsf,
					   double *ldirty, double *psf,
					   const int nscales, const int nmoments,
					   const int ny, const int nx,
					   double *windowstack,
					   double gain, double absolutethresh, int niter,
					   const char* findpeak = "CASA");

int msmfsclean_kernel_c(PyArrayObject *&m_model_obj,
                 PyArrayObject *&residual_obj,
                 PyArrayObject *&scalestack_obj,
                 PyArrayObject *&smresidual_obj,
                 PyArrayObject *&ssmmpsf_obj,
                 PyArrayObject *&hsmmpsf_obj,
                 PyArrayObject *&ihsmmpsf_obj,
                 PyArrayObject *&ldirty_obj,
                 PyArrayObject *&psf_obj,
                 PyArrayObject *&windowstack_obj,
                 double gain,
                 double absolutethresh,
                 int niter,
                 const char* findpeak);

#endif