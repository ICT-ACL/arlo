#include <map>
#include <set>
#include <stdlib.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <algorithm>
#include <string>
#include <omp.h>
// #include <cblas.h>

// #include "fftw.h"
#include "cleaners.h"

using namespace std;


// /**
// * Evaluates the PROLATE SPHEROIDAL WAVEFUNCTION
// * m=6, alpha = 1 from Schwab, Indirect Imaging (1984).
// * This is one factor in the basis function.
// * 
// * Code adapted Anna's f90 PROFILE (gridder.f90) code
// * which was adapted from Tim Cornwell's C++ SphFuncVisGridder
// * developed for CONRAD for ASKAP. **This seems to be commented
// * out of the currect ASKAPsoft code... not sure why**
// * Stole this back from Anna!
// */
// double spheroidal_function(double vnu) {
//  int n_p = 4, n_q = 2;

    
//  int p[] = {8.203343e-2, -3.644705e-1, 6.278660e-1, -5.335581e-1, 2.312756e-1, 
//             4.028559e-3, -3.697768e-2, 1.021332e-1, -1.201436e-1, 6.412774e-2};

//  int q[] = {1.0000000, 8.212018e-1, 2.078043e-1, 
//             1.0000000, 9.599102e-1, 2.918724e-1};

//  double value = 0.0;

//  int part, nuend;
//  if (vnu >= 0.0 && vnu < 0.75) {
//      part = 0;
//      nuend = 0.75;
//  }
//  else if (vnu >= 0.75 && vnu <= 1.0) {
//      part = 1;
//      nuend = 1.0;
//  }
//  else {
//      value = 0.0;
//      return value;
//  }

//  int top = p[part*5], bot = q[part*3];
//  double delnusq = vnu * vnu - nuend * nuend;

//  for (int k = 1; k <= n_p; k++) {
//      double factor = pow(delnusq, k);
//      top += p[part*5+k] * factor;
//  }

//  for (int k = 1; k <= n_q; k++) {
//      double factor = pow(delnusq, k);
//      bot += q[part*3+k] * factor;
//  }

//  if (bot != 0.0) {
//      value = top / bot;
//  }
//  else {
//      value = 0.0;
//  }

//  if (value < 0.0) {
//      value = 0.0;
//  }

//  return value;
// }
    
    
// /**
// *
// * Create a cube consisting of the scales
// * param scale_shape: desired shape of stack
// * param scales: scales (in pixels)
// * param norm: Normalise each plane to unity?
// * return: stack
// */
// void create_scalestack(double *basis, const int*scale_shape, int *scales, bool norm) {
    
//  int nscales = scale_shape[0], nx = scale_shape[1], ny = scale_shape[2];
    
//  int xcen = int(ceil(double(nx) / 2.0));
//  int ycen = int(ceil(double(ny) / 2.0));

//  for (int iscale = 0; iscale < nscales; i++) {
//      double *basis_local = basis + iscale * nx * ny;

//      int halfscale = int(ceil(scales[iscale] / 2.0));
//      if (scales[iscale] > 0.0):
//          double rscale = 1.0 / (double(scales[iscale]) / 2.0);
//          double rscale2 = rscale * rscale;

//          // int begin = xcen - halfscale - 1, end = xcen + halfscale + 1;
//          // int len = end - begin;
//          // double *xs = (double*)malloc(sizeof(double) * len);
//          // double *fx = (double*)malloc(sizeof(double) * len);

//          // for (int i = 0; i < len; i++) {
//          //  xs[i] = begin + i;
//          //  fxs[i] = begin + i - xcen;
//          // }

//          for (int y = ycen - halfscale - 1; y < ycen + halfscale + 1; y++) {
//              for (int x = xcen - halfscale - 1; x < xcen + halfscale + 1; x++) {
//                  double fx = x - xcen, fy = y - ycen;
//                  double r2 = rscale2 * (fx * fx + fy * fy);
//                  double r = sqrt(r2);
//                  double val = spheroidal_function(r) * (1.0 - r * r);
//                  basis_local[x * ny + y] = min(val, 0.0);
//              }
//          }

//          if (norm) {
//              double sum = 0.0;
//              for (int i = 0; i < nx * ny; i++) {
//                  sum += basis_local[i];
//              }

//              for (int i = 0; i < nx * ny; i++) {
//                  basis_local[i] /= sum;
//              }
//          }

//      else {
//          basis[xcen * ny + ycen] = 1.0;
//      }
//  }
// }


// /*
//     """Convolve img by the specified scalestack, returning the resulting stack

//     :param scalestack: stack containing the scales
//     :param img: Image to be convolved
//     :return: stack
//     """
// */
// void convolve_scalestack(double *convolved, double *scalestack, double *img, const int *scale_shape) {
//  int nscales = scale_shape[0], nx = scale_shape[1], ny = scale_shape[2];

//  complex_t *tmp_img = (complex_t*)malloc(sizeof(complex_t) * nx * ny);
//  complex_t *ximg = (complex_t*)malloc(sizeof(complex_t) * nx * ny);

//  fftshift(tmp_img, img, nx, ny);
//  fftw_fft2(tmp_img, tmp_img, nx, ny);
//  fftshift(ximg, tmp_img, nx, ny);

//  complex_t *tmp_scale = (complex_t*)malloc(sizeof(complex_t) * nx * ny);
//  complex_t *xscale = (complex_t*)malloc(sizeof(complex_t) * nx * ny);
//  complex_t *xmult = (complex_t*)malloc(sizeof(complex_t) * nx * ny);
//  complex_t *tmp_xmult = (complex_t*)malloc(sizeof(complex_t) * nx * ny);
//  complex_t *tmp_convolved = (complex_t*)malloc(sizeof(complex_t) * nx * ny);

//  for (int iscale = 0; iscale < nscales; i++) {
//      fftshift(tmp_scale, scalestack + iscale * nx * nx, ny, ny);
//      fftw_fft2(tmp_scale, tmp_scale, nx, ny);
//      fftshift(xscale, tmp_scale, nx, xy);

//      for (int i = 0; i < nx * ny; i++) {
//          xmult[i] = ximg[i] * conj(xscale);
//      }

//      ifftshift(tmp_xmult, xmult, nx, ny);
//      fftw_ifft2(tmp_xmult, tmp_xmult, nx, ny);
//      ifftshift(tmp_convolved, tmp_xmult, nx, ny);

//      double *convolved_local = convolved + iscale * nx * ny;
//      for (int i = 0; i < nx * ny; i++) {
//          convolved_local[i] = creal(tmp_convolved[i]);
//      }
//  }

//  free(tmp_img);
//  free(ximg);
//  free(tmp_scale);
//  free(xscale);
//  free(xmult);
//  free(tmp_xmult);
//  free(tmp_convolved);
// }


// /**
// *
// * Calculate scale-dependent moment residuals
// * Part of the initialisation for Algorithm 1: lines 12 - 17
// * :param residual: residual [nmoments, nx, ny]
// * :return: scale-dependent moment residual [nscales, nmoments, nx, ny]
// */
// void calculate_scale_moment_residual(double *scale_moment_residual,
//                                      double *residual, const int *residual_shape,
//                                      double *scalestack, const int *scale_shape) {
//  int nmoments = residual_shape[0], nx = residual_shape[1], ny = residual_shape[2];
//  int nscales = scale_shape[0];

//  double *convolved = (double*)malloc(sizeof(double) * nscales * nx * ny);

//  for (int t = 0; t < nmoments; t++) {
//      convolve_scalestack(convolved, scalestack, residual + t * nx * ny, scale_shape);
//      for (int iscale = 0; iscale < nscales; i++) {
//          memcpy(scale_moment_residual + (iscale * nmoments + t) * nx * ny,
//                 convolved + iscale * nx * ny,
//                 sizeof(double) * nx * ny);
//      }
//  }

//  free(convolved);
// }


// /**
// """Convolve img by the specified scalestack, returning the resulting stack

//     :param scalestack: stack containing the scales
//     :param img: Image to be convolved
//     :return: Twice convolved image [nscales, nscales, nx, ny]
//     """
// */
// void convolve_convolve_scalestack(double *convolved, double *scalestack, double *img, const int *scale_shape) {
//  int nscales = scale_shape[0], nx = scale_shape[1], ny = scale_shape[2];
//  int convolved_shape = nscales * nscales * nx * ny;
//  int xscale_shape = nscales * nx * ny;

//  complex_t *tmp_img = (complex_t*)malloc(sizeof(complex_t) * nx * ny);
//  complex_t *ximg = (complex_t*)malloc(sizeof(complex_t) * nx * ny);

//  fftshift(tmp_img, img, nx, ny);
//  fftw_fft2(tmp_img, tmp_img, nx, ny);
//  fftshift(ximg, tmp_img, nx, ny);

//  complex_t *xscale = (complex_t*)malloc(sizeof(complex_t) * xscale_shape);
//  complex_t *tmp_stack = (complex_t*)malloc(sizeof(complex_t) * nx * ny);

//  for (int s = 0; s < nscales; s++) {
//      fftshift(tmp_stack, scalestack + s * nx * nx, ny, ny);
//      fftw_fft2(tmp_stack, tmp_stack, nx, ny);
//      fftshift(xscale + s * nx * ny, tmp_stack, nx, ny);
//  }

//  complex_t *xmult = (complex_t*)malloc(sizeof(complex_t) * nx * ny);
//  complex_t *tmp_xmult = (complex_t*)malloc(sizeof(complex_t) * nx * ny);
//  complex_t *tmp_convolved = (complex_t*)malloc(sizeof(complex_t) * nx * ny);

//  for (int s = 0; s < nscales; s++) {
//      for (int p = 0; p < nscales; p++) {

//          for (int i = 0; i < nx * ny; i++) {
//              xmult[i] = ximg[i] * xscale[p*nx*ny + i] * conj(xscale[s*nx*ny + i]);
//          }

//          ifftshift(tmp_xmult, xmult, nx, ny);
//          fftw_ifft2(tmp_xmult, tmp_xmult, nx, ny);
//          ifftshift(tmp_convolved, tmp_xmult, nx, ny);

//          double *convolved_local = convolved + (s * nscales + p) * nx * ny;
//          for (int i = 0; i < nx * ny; i++) {
//              convolved_local[i] = creal(tmp_convolved[i]);
//          }
//      }
//  }

//  free(tmp_img);
//  free(ximg);
//  free(xscale);
//  free(tmp_stack)
//  free(xmult);
//  free(tmp_xmult);
//  free(tmp_convolved);
// }


// /**
// * Calculate scale-dependent moment psfs
// * Part of the initialisation for Algorithm 1
// * :param psf: psf
// * :return: scale-dependent moment psf [nscales, nscales, nmoments, nmoments, nx, ny]
// */
// void calculate_scale_scale_moment_moment_psf(double *scale_scale_moment_moment_psf,
//                                          double *psf, double *scalestack, 
//                                          const int *psf_shape, const int *scale_shape) {
    
//  int nmoments2 = psf_shape[0], nx = psf_shape[1], ny = psf_shape[2];
//  int nmoments = nmoments2 / 2;
//  int nscales = scale_shape[0], snx = scale_shape[1], sny = scale_shape[2];

//  double *convolved = (double*)malloc(sizeof(double) * nscales * nscales * nx * ny);

//  for (int t = 0; t < nmoments; t++) {
//      for (int q = 0; q < nmoments; q++) {

//          convolve_convolve_scalestack(convolved, scalestack, 
//                                       psf + (t+q) * nx * ny, scale_shape);
            
//          for (int i = 0; i < nscales; i++) {
//              for (int j = 0; j < nscales; j++) {
//                  memcpy(scale_scale_moment_moment_psf + (((i*nscales+j) * nmoments + t) * nmoments + q),
//                      convolved + i * nscales + j,
//                      sizeof(double) * nx * ny);

//              }
//          }
//      }
//  }

//  free(convolved);
// }


// void calculate_scale_inverse_moment_moment_hessian(double *scale_moment_moment_hessian, 
//                                                 double *scale_inverse_moment_moment_hessian, 
//                                                 double *scale_scale_moment_moment_psf,
//                                                 const int *ssmmpsf_shape) {
//  int nscales = ssmmpsf_shape[0], nmoments = ssmmpsf_shape[2], nx = ssmmpsf_shape[4], ny = ssmmpsf_shape[5];

//  for (int s = 0; s < nscales; s++) {
//      for (int i = 0; i < nmoments * nmoments; i++) {
//          scale_moment_moment_hessian[s * nmoments * nmoments + i] = 
//              scale_scale_moment_moment_psf[((((s*nscales)+s)*nmoments*nmoments+i)*nx+nx/2)*ny+ny/2];
//      }
//  }
// }



// /**
//     """ Perform image plane multiscale multi frequency clean

//     This algorithm is documented as Algorithm 1 in: U. Rau and T. J. Cornwell, “A multi-scale multi-frequency
//     deconvolution algorithm for synthesis imaging in radio interferometry,” A&A 532, A71 (2011). Note that
//     this is only the image plane parts.

//     Specific code is linked to specific lines in that algorithm description.

//     This version operates on numpy arrays that have been converted to moments on the last axis.

//     :param fracthresh:
//     :param dirty: The dirty image, i.e., the image to be deconvolved
//     :param psf: The point spread-function
//     :param window: Regions where clean components are allowed. If True, all of the dirty image is allowed
//     :param gain: The "loop gain", i.e., the fraction of the brightest pixel that is removed in each iteration
//     :param thresh: Cleaning stops when the maximum of the absolute deviation of the residual is less than this value
//     :param niter: Maximum number of components to make if the threshold "thresh" is not hit
//     :param scales: Scales (in pixels width) to be used
//     :param fracthres: Fractional stopping threshold
//     :param ntaylor: Number of Taylor terms
//     :param findpeak: Method of finding peak in mfsclean: 'Algorithm1'|'CASA'|'ARL', Default is ARL.
//     :return: clean component image, residual image
//     """
// */    
// void msmfsclean(double *m_model, double *residual,
//              double *dirty, const int *dirty_shape,
//              double *psf, const int *psf_shape,
//              double *window, const int *window_shape,
//              float gain, float thresh, int niter,
//              int *scales, const int nscales,
//              float fracthresh, string findpeak) {
    
//  int nmoments_d = dirty_shape[0], nx_d = dirty_shape[1], ny_d =dirty_shape[2];
//  int nmoments_p = psf_shape[0], nx_p = psf_shape[1], ny_p = psf_shape[2];

//  int dirty_len = nmoments_d * nx_d * ny_d;
//  int psf_len = nmoments_p * nx_p * ny_p;

//  double pmax = psf[0];
//  for (int i = 1; i < psf_len; i++) {
//      pmax = pmax > psf[i] ? pmax : psf[i];
//  }

//  double *lpsf = (double*)malloc(sizeof(double) * psf_len);
//  double *ldirty = (double*)malloc(sizeof(double) * dirty_len);

//  for (int i = 0; i < psf_len; i++) {
//      lpsf[i] = psf[i] / pmax;
//  }

//  for (int i = 0; i < dirty_len; i++) {
//      ldirty[i] = dirty[i] / pmax;
//  }

//  int scale_shape[] = {nscales, nx_d, ny_d};

//  double *scalestack = (double*)malloc(sizeof(double) * nscales * nx_d * ny_d);
//  create_scalestack(scalestack, scale_shape, scales, true);

//  double *smresidual = (double*)malloc(sizeof(double) * nscales * nmoments_p * nx_p * ny_p);
//  calculate_scale_moment_residual(smresidual, ldirty, dirty_shape, scalestack, scale_shape);

//  double *ssmmpsf = (double*)malloc(sizeof(double) * nscales * nscales * nmoments_p * nmoments_p * nx_p * ny_p);
//  calculate_scale_scale_moment_moment_psf(ssmmpsf, lpsf, scalestack, psf_shape, scale_shape);

//  double *hsmmpsf = (double*)malloc(sizeof(double) * nscales * nmoments_p, nmoments_p);
//  double * ihsmmpsf = (double*)malloc(sizeof(double) * nscales * nmoments_p, nmoments_p);
//  calculate_scale_inverse_moment_moment_hessian(ssmmpsf, nscales, nmoments_p);



//  free(lpsf);
//  free(ldirty);
//  free(scalestack);
//  free(smresidual);
//  free(ssmmpsf);
//  free(hsmmpsf);
//  free(ihsmmpsf);
// }



/**
   """ Find the indices where two arrays overlap

    :param a1: First array
    :param a2: Second array
    :param shiftx: Shift in x applied to a1
    :param shifty: Shift in y applied to a2
    :return: (limits in a1, limits in a2)
    """
*/    
void overlap_indices(int *lhs, int *rhs, 
                    double *res, double *psf, 
                    int peakx, int peaky,
                    const int *res_shape, const int *psf_shape) {

    int nx = res_shape[0], ny = res_shape[1];
    int psfwidthx = psf_shape[0] / 2, psfwidthy = psf_shape[1] / 2;
    int psfpeakx = psf_shape[0] / 2, psfpeaky = psf_shape[1] / 2;

    int res_lower[] = {max(0, peakx - psfwidthx), max(0, peaky - psfwidthy)};
    int res_upper[] = {min(nx, peakx + psfwidthx), min(peaky + psfwidthy, ny)};

    int psf_lower[] = {max(0, psfpeakx + (res_lower[0] - peakx)), max(0, psfpeaky + (res_lower[1] - peaky))};
    int psf_upper[] = {min(psf_shape[0], psfpeakx + (res_upper[0] - peakx)), min(psfpeaky + (res_upper[1] - peaky), psf_shape[1])};

    lhs[0] = res_lower[0], lhs[1] = res_upper[0], lhs[2] = res_lower[1], lhs[3] = res_upper[1];
    rhs[0] = psf_lower[0], rhs[1] = psf_upper[0], rhs[2] = psf_lower[1], rhs[3] = psf_upper[1];
}



/*
    """ Calculate the principal solution in moment space for each scale

    Lines 20 - 26

    :param smresidual: scale-dependent moment residual [nscales, nmoments, nx, ny]
    :param imhsmmpsf: Inverse of scale dependent moment moment Hessian
    :return: Decoupled residual images [nscales, nmoments, nx, ny]
    """
    # ihsmmpsf: nscales, nmoments, nmoments
    # smresidual: nscales, nmoments, nx, ny
*/
void calculate_scale_moment_principal_solution(double *smpsol, double *smresidual, double *ihsmmpsf,
                                               const int nmoments, const int nscales, 
                                               const int nx, const int ny) {

    const int nz = nx * ny;
    memset(smpsol, 0, sizeof(double) * nscales * nmoments * nz);

    // for (int s = 0; s < nscales; s++) {
    //     for (int n = 0; n < nmoments; n++) {
    //         for (int x = 0; x < nx; x++) {
    //             for (int y = 0; y < ny; y++) {
 //                 for (int m = 0; m < nmoments; m++) {
    //                     smpsol[(((s*nmoments+n)*nx+x)*ny+y)] += ihsmmpsf[((s*nmoments+m)*nmoments+n)] * 
    //                                                             smresidual[(((s*nmoments+m)*nx+x)*ny+y)];
    //                 }
    //             }
    //         }
    //     }
    // }

    #pragma omp parallel
    {
        int nthreads = omp_get_max_threads();
        int tid = omp_get_thread_num();
        int low = tid * nz / nthreads;
        int high = (tid + 1) * nz / nthreads;

        for (int s = 0; s < nscales; s++) {
            for (int n = 0; n < nmoments; n++) {
                for (int m = 0; m < nmoments; m++) {

                    double v = ihsmmpsf[(s*nmoments+m)*nmoments+n];
                    double *smpsol_m = smpsol + (s * nmoments + n) * nz;
                    double *smresidual_m = smresidual + (s * nmoments + m) * nz;

                    # pragma simd
                    for (int z = low; z < high; z++) {
                        smpsol_m[z] += v * smresidual_m[z];
                    }

                }
            }
        }
    }
}



/*
    """Find the optimum scale for moment zero

    Line 27 of Algorithm 1

    :param smpsol: Decoupled residual images for each scale and moment
    :return: x, y, optimum scale for peak
    """
*/    
void find_optimum_scale_zero_moment(int &mx, int &my, int &mscale, double *smpsol, double *windowstack, 
                                    const int nmoments, const int nscales, const int nx, const int ny) {

    double optimum = 0.0;
    int index = 0;

    if (windowstack != NULL) {
        for (int s = 0; s < nscales; s++) {
            int begin_s = s * nmoments * nx * ny;
            int begin_w = s * nx * ny;

            for (int z = 0; z < nx * ny; z++) {
                double val = smpsol[begin_s + z] * windowstack[begin_w + z];
                if (val > optimum) {
                    optimum = val;
                    index = begin_s + z;
                }
            }
        }
    }
    else {
        for (int s = 0; s < nscales; s++) {
            int begin_s = s * nmoments * nx * ny;

            for (int z = 0; z < nx * ny; z++) {
                double val = smpsol[begin_s + z];
                if (val > optimum) {
                    optimum = val;
                    index = begin_s + z;
                }
            }
        }
    }

    mscale = index / (nmoments*nx*ny);
    index %= (nmoments * nx * ny);
    mx = index / ny;
    my = index % ny;
}




/*
    """Find the optimum peak using one of a number of algorithms

    """
*/    
void find_global_optimum(int &mscale, int &mx, int &my, double *mval, 
                        double *hsmmpsf, double *ihsmmpsf, 
                        double *smresidual, double *windowstack, 
                        const int nmoments, const int nscales, 
                        const int nx, const int ny,
                        const char *findpeak,
                        double &time1, double &time2, double &time3) {

    double time;
    double *smpsol = (double*)malloc(sizeof(double) * nscales * nmoments * nx * ny);
    time = -omp_get_wtime();
    calculate_scale_moment_principal_solution(smpsol, smresidual, ihsmmpsf, nmoments, nscales, nx, ny);
    time += omp_get_wtime();
    time1 += time;

    if (strcmp(findpeak, "Algorithm1") == 0) {
        find_optimum_scale_zero_moment(mx, my, mscale, smpsol, windowstack, nmoments, nscales, nx, ny);
    }

    else if (strcmp(findpeak, "CASA") == 0) {
        double *dchisq = (double*)malloc(sizeof(double) * nscales * 1 * nx * ny);

        for (int scale = 0; scale < nscales; scale++) {

            double *dchisq_s = dchisq + scale * 1 * nx * ny;
            double *smpsol_s = smpsol + scale * nmoments * nx * ny;
            double *smresidual_s = smresidual + scale * nmoments * nx * ny;
            double *hsmmpsf_s = hsmmpsf + scale * nmoments * nmoments;

            for (int monent1 = 0; monent1 < nmoments; monent1++) {

                double *smpsol_sm = smpsol_s + monent1 * nx * ny;
                double *smresidual_sm = smresidual_s + monent1 * nx * ny;

                for (int i = 0; i < nx; i++) {
                    for (int j = 0; j < ny; j++) {
                        dchisq_s[i * ny + j] += 2.0 * smpsol_sm[i * ny + j] * smresidual_sm[i * ny + j];
                    }
                }

                for (int monent2 = 0; monent2 < nmoments; monent2++) {

                    for (int i = 0; i < nx; i++) {
                        for (int j = 0; j < ny; j++) {
                            dchisq_s[i * ny + j] -= hsmmpsf_s[monent1 * nmoments + monent2] * 
                                                    smpsol_sm[i * ny + j] * 
                                                    smpsol_s[(monent2 * nx + i) * ny + j];
                        }
                    }
                }
            }
        }

        find_optimum_scale_zero_moment(mx, my, mscale, dchisq, windowstack, nmoments, nscales, nx, ny);
        
        free(dchisq);
    }

    else {
        double *prod_sol_res = (double*)malloc(sizeof(double) * nscales * nmoments * nx * ny);

        time = -omp_get_wtime();
        # pragma omp parallel for
        for (int i = 0; i < nscales * nmoments * nx * ny; i++) {
            prod_sol_res[i] = smpsol[i] * smresidual[i];
        }
        time += omp_get_wtime();
        time2 += time;

        time = -omp_get_wtime();
        find_optimum_scale_zero_moment(mx, my, mscale, prod_sol_res, windowstack, nmoments, nscales, nx, ny);
        time += omp_get_wtime();
        time3 += time;

        free(prod_sol_res);
    }

    double *smpsol_local = smpsol + mscale * nmoments * nx * ny;
    for (int i = 0; i < nmoments; i++) {
        mval[i] = smpsol_local[(i * nx + mx) * ny + my];
    }

    free(smpsol);
}




/**
    """Update model with an appropriately scaled and centered blob for each moment

    """
*/    
void update_moment_model(double *m_model, double *scalestack,
                        int* lhs, int *rhs, 
                        double gain, int mscale, double *mval,
                        const int nmoments, const int nscales, 
                        const int nx, const int ny) {

    int ldx = lhs[1] - lhs[0], ldy = lhs[3] - lhs[2];
    int rdx = rhs[1] - rhs[0], rdy = rhs[3] - rhs[2];
    assert(ldx == rdx && ldy == rdy);
    int dx = ldx, dy = ldy;

    for (int t = 0; t < nmoments; t++) {

        double val = gain * mval[t];
        double *m_model_local = m_model + t * nx * ny ;
        double *scalestack_local = scalestack + mscale * nx * ny;

        for (int i = 0; i < dx; i++) {
            for (int j = 0; j < dy; j++) {

                m_model_local[(lhs[0]+i) * ny + (lhs[2]+j)] += scalestack_local[(rhs[0]+i) * ny + (rhs[2]+j)] * val;
            
            }
        }
    }
}



/*
    """ Update residual by subtracting the effect of model update for each moment

    """
*/    
void update_scale_moment_residual(double *smresidual, double *ssmmpsf, 
                                int *lhs, int *rhs, 
                                double gain, int mscale, double *mval, 
                                const int nmoments, const int nscales, 
                                const int nx, const int ny) {

    int ldx = lhs[1] - lhs[0], ldy = lhs[3] - lhs[2];
    int rdx = rhs[1] - rhs[0], rdy = rhs[3] - rhs[2];
    assert(ldx == rdx && ldy == rdy);
    int dx = ldx, dy = ldy;
    
    double *ssmmpsf_0 = ssmmpsf + mscale * nscales * nmoments * nmoments * nx * ny;

    // for (int mi = 0; mi < nscales; mi++) {
    //  for (int mj = 0; mj < nmoments; mj++) {

    //      double *ssmmpsf_local = ssmmpsf_m + (mi * nmoments + mj) * nmoments * nx * ny;
    //      double *smresidual_local = smresidual + (mi * nmoments + mj) * nx * ny;

    //      for (int i = 0; i < dx; i++) {
    //          for (int j = 0; j < dy; j++) {

    //              double sum = 0.0;
    //              for (int v = 0; v < nmoments; v++) {
    //                  sum += ssmmpsf_local[(v * nx + (rhs[0]+i)) * ny + (rhs[2]+j)] * mval[v];
    //              }

    //              smresidual_local[(lhs[0]+i) * ny + (lhs[2]+j)] -= gain * sum;

    //          }
    //      }
    //  }
    // }

    #pragma omp parallel
    {
        int nthreads = omp_get_max_threads();
        int tid = omp_get_thread_num();
        int low = tid * dx / nthreads;
        int high = (tid + 1) * dx / nthreads;

        for (int s = 0; s < nscales; s++) {
            for (int n = 0; n < nmoments; n++) {

                double *ssmmpsf_n = ssmmpsf_0 + (s * nmoments + n) * nmoments * nx * ny;
                double *smresidual_n = smresidual + (s * nmoments + n) * nx * ny;

                for (int x = low; x < high; x++) {

                    #pragma simd
                    for (int y = 0; y < dy; y++) {

                        double sum = 0.0;
                        for (int m = 0; m < nmoments; m++) {
                            sum += ssmmpsf_n[(m * nx + (rhs[0]+x)) * ny + (rhs[2]+y)] * mval[m];
                        }

                        smresidual_n[(lhs[0]+x) * ny + (lhs[2]+y)] -= gain * sum;

                    }
                }
            }
        }
    }
}


/**
    """ Perform image plane multiscale multi frequency clean

    This algorithm is documented as Algorithm 1 in: U. Rau and T. J. Cornwell, “A multi-scale multi-frequency
    deconvolution algorithm for synthesis imaging in radio interferometry,” A&A 532, A71 (2011). Note that
    this is only the image plane parts.

    Specific code is linked to specific lines in that algorithm description.

    This version operates on numpy arrays that have been converted to moments on the last axis.

    :param fracthresh:
    :param dirty: The dirty image, i.e., the image to be deconvolved
    :param psf: The point spread-function
    :param window: Regions where clean components are allowed. If True, all of the dirty image is allowed
    :param gain: The "loop gain", i.e., the fraction of the brightest pixel that is removed in each iteration
    :param thresh: Cleaning stops when the maximum of the absolute deviation of the residual is less than this value
    :param niter: Maximum number of components to make if the threshold "thresh" is not hit
    :param scales: Scales (in pixels width) to be used
    :param fracthres: Fractional stopping threshold
    :param ntaylor: Number of Taylor terms
    :param findpeak: Method of finding peak in mfsclean: 'Algorithm1'|'CASA'|'ARL', Default is ARL.
    :return: clean component image, residual image
    """
*/    
void msmfsclean_kernel(double *m_model, double *residual,
                       double *scalestack, double *smresidual,
                       double *ssmmpsf, double *hsmmpsf, double *ihsmmpsf,
                       double *ldirty, double *psf,
                       const int nscales, const int nmoments,
                       const int nx, const int ny,
                       double *windowstack,
                       double gain, double absolutethresh, int niter,
                       const char *findpeak) {

    // m_model: nmoments, nx, ny
    // dirty: nmoments, nx, ny
    // ldirty: nmoments, nx, ny
    // psf: 2*nmoments, nx, ny
    // lpsf: 2*nmoments, nx, ny
    // smresidual: nscales, nmoments, nx, ny
    // scalestack: nscales, nx, ny
    // pscalestack: nscales, nx, ny
    // ssmmpsf: nscales, nscales, nmoments, nmoments, nx, ny
    // hsmmpsf: nscales, nmoments, nmoments
    // ihsmmpsf: nscales, nmoments, nmoments
    const int dirty_len = nmoments * nx * ny;
    const int psf_len = 2 * nmoments * nx * ny;
    const int scale_len = nscales * nx * ny;
    const int hsmm_len = nscales * nmoments * nmoments;
    const int ssmm_len = nscales * nscales * nmoments * nmoments * nx * ny;
    const int res_len = nscales * nmoments * nx * ny;

    double pmax = psf[0];
    for (int i = 1; i < psf_len; i++) {
        pmax = max(pmax, psf[i]);
    }

    int mscale, mx, my;
    double *mval = (double*)malloc(sizeof(double) * nmoments);
    int lhs[4], rhs[4];

    int ldirty0_shape[] = {nx, ny};
    int psf0_shape[] = {nx, ny};

    double time1 = 0, time2 = 0, time3 = 0, time4 = 0, time;
    double timea = 0, timeb = 0, timec = 0;
    for (int i = 0; i < niter; i++) {
#ifdef DEBUG
        if (i % 100 == 0) {
            printf("iter = %d\n", i);
        }
#endif
        time = -omp_get_wtime();
        find_global_optimum(mscale, mx, my, mval, 
                            hsmmpsf, ihsmmpsf, smresidual, windowstack, 
                            nmoments, nscales, nx, ny, findpeak,
                            timea, timeb, timec);
        time += omp_get_wtime();
        time1 += time;

        double mval_max = mval[0];
        for (int v = 1; v < nmoments; v++) {
            mval_max = max(mval_max, mval[v]);
        }
        if (fabs(mval_max) < absolutethresh) {
            break;
        }

        overlap_indices(lhs, rhs, ldirty, psf, mx, my, 
                        ldirty0_shape, psf0_shape);

        time = -omp_get_wtime();
        update_moment_model(m_model, scalestack, lhs, rhs, gain, mscale, mval, 
                            nmoments, nscales, nx, ny);
        time += omp_get_wtime();
        time3 += time;

        time = -omp_get_wtime();
        update_scale_moment_residual(smresidual, ssmmpsf, lhs, rhs, gain, mscale, mval, 
                                     nmoments, nscales, nx, ny);
        time += omp_get_wtime();
        time4 += time;
    }

    # pragma omp parallel for
    for (int i = 0; i < nmoments * nx * ny; i++) {
        residual[i] = pmax * smresidual[i];
    }

    printf("time find_global_optimum   : %.4lf\n", time1);
    printf("-- time calculate_solution : %.4lf\n", timea);
    printf("-- time smpsol * smresidual: %.4lf\n", timeb);
    printf("-- time find_zero_moment   : %.4lf\n", timec);
    printf("time update_model          : %.4lf\n", time3);
    printf("time update_residual       : %.4lf\n", time4);

    free(mval);
}





// void msmfsclean_kernel(double *m_model, double *residual,
//                        double *scalestack, double *smresidual,
//                        double *ssmmpsf, double *hsmmpsf, double *ihsmmpsf,
//                        double *ldirty, double *psf,
//                        const int nscales, const int nmoments,
//                        const int ny, const int nx,
//                        double *windowstack,
//                        double gain, double absolutethresh, int niter,
//                        const char* findpeak) {
//     printf("success\n");
// }


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
                 const char* findpeak) {

    int nscales = smresidual_obj->dimensions[0];
    int nmoments = smresidual_obj->dimensions[1];
    int ny = smresidual_obj->dimensions[2];
    int nx = smresidual_obj->dimensions[3];

    double *m_model = (double*)m_model_obj->data;
    double *residual = (double*)residual_obj->data;
    double *scalestack = (double*)scalestack_obj->data;
    double *smresidual = (double*)smresidual_obj->data;
    double *ssmmpsf = (double*)ssmmpsf_obj->data;
    double *hsmmpsf = (double*)hsmmpsf_obj->data;
    double *ihsmmpsf = (double*)ihsmmpsf_obj->data;
    double *ldirty = (double*)ldirty_obj->data;
    double *psf = (double*)psf_obj->data;
    double *windowstack_tmp = (double*)windowstack_obj->data;
    double *windowstack;

    double sum = 0.0;
    for (int i = 0; i < nscales * nx * ny; i++) {
        sum += windowstack[i];
    }
    if (fabs(sum) < 1e-6) {
        windowstack = NULL;
    }
    else {
        windowstack = windowstack_tmp;
    }

    msmfsclean_kernel(m_model, residual,
                       scalestack, smresidual,
                       ssmmpsf, hsmmpsf, ihsmmpsf,
                       ldirty, psf,
                       nscales, nmoments,
                       ny, nx,
                       windowstack,
                       gain, absolutethresh, niter,
                       findpeak);

    return 0;
}