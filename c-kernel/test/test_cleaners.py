import os
import sys
sys.path.append('../../arl-python')

import numpy as np
import time

from arl.image.cleaners import *
from utils import *


def msmfsclean_simplify(dirty, psf, window, gain, thresh, niter, scales, fracthresh, findpeak='CASA'):

    assert 0.0 < gain < 2.0
    assert niter > 0
    assert len(scales) > 0

    m_model = np.zeros(dirty.shape)

    nscales = len(scales)

    pmax = psf.max()
    assert pmax > 0.0

    psfpeak = np.argmax(np.fabs(psf))
    dmax = dirty.max()
    dpeak = np.argmax(dirty)
    lpsf = psf / pmax
    ldirty = dirty / pmax

    nmoments, ny, nx = dirty.shape
    assert psf.shape[0] == 2 * nmoments

    # Create the "scale basis functions" in Algorithm 1
    scaleshape = [nscales, ldirty.shape[1], ldirty.shape[2]]
    scalestack = create_scalestack(scaleshape, scales, norm=True)

    pscaleshape = [nscales, lpsf.shape[1], lpsf.shape[2]]
    pscalestack = create_scalestack(pscaleshape, scales, norm=True)

    # Calculate scale convolutions of moment residuals
    smresidual = calculate_scale_moment_residual(ldirty, scalestack)
    smresidual0 = smresidual.copy()

    # Calculate scale scale moment moment psf, Hessian, and inverse of Hessian
    # scale scale moment moment psf is needed for update of scale-moment residuals
    # Hessian is needed in calculation of optimum for any iteration
    # Inverse Hessian is needed to calculate principal solution in moment-space
    ssmmpsf = calculate_scale_scale_moment_moment_psf(lpsf, pscalestack)
    hsmmpsf, ihsmmpsf = calculate_scale_inverse_moment_moment_hessian(ssmmpsf)

    # The window is scale dependent - we form it by smoothing and thresholding
    # the input window. This prevents components being placed too close to the
    # edge of the Image.

    if window is None:
        windowstack = None
    else:
        windowstack = np.zeros_like(scalestack)
        windowstack[convolve_scalestack(scalestack, window) > 0.9] = 1.0

    absolutethresh = max(thresh, fracthresh * np.fabs(smresidual[0, 0, :, :]).max())

    # Start iterations
    scale_counts = np.zeros(nscales, dtype='int')
    scale_flux = np.zeros(nscales)

    # Use original algorithm
    start = time.time()

    for i in range(niter):
        # Find the optimum scale and location.
        mscale, mx, my, mval = find_global_optimum(hsmmpsf, ihsmmpsf, smresidual, windowstack, findpeak)

        scale_counts[mscale] += 1
        scale_flux[mscale] += mval[0]

        # Are we ready to stop yet?
        peak = np.max(np.fabs(mval))
        if peak < absolutethresh:
            break

        # Calculate indices needed for lhs and rhs of updates to model and residual
        lhs, rhs = overlapIndices(ldirty[0, ...], psf[0, ...], mx, my)

        m_model = update_moment_model(m_model, pscalestack, lhs, rhs, gain, mscale, mval)

        smresidual = update_scale_moment_residual(smresidual, ssmmpsf, lhs, rhs, gain, mscale, mval)

    residual = pmax * smresidual[0, :, :, :]

    stop = time.time()
    print('Original Time:  {:.2f}s'.format(stop - start))

    return m_model, residual, pscalestack, smresidual0, \
            ssmmpsf, hsmmpsf, ihsmmpsf, ldirty, psf



def test_cleaners(data_dir):

    nscales, nmoments, nx, ny = 4, 3, 512, 512
    dirty = create_random_data((nmoments, ny, nx), -1000, 1000, 'float')
    psf = create_random_data((nmoments*2, ny, nx), -5, 5, 'float')

    m_model, residual, pscalestack, smresidual0, \
    ssmmpsf, hsmmpsf, ihsmmpsf, ldirty, psf \
        = msmfsclean_simplify(dirty, psf, None, gain=0.7, thresh=0.01, niter=2, scales=[0, 3, 10, 30],\
                        fracthresh=0.001, findpeak='ARL')


    store_data(os.path.join(data_dir, 'm_model.dat'), m_model)
    store_data(os.path.join(data_dir, 'residual.dat'), residual)
    store_data(os.path.join(data_dir, 'pscalestack.dat'), pscalestack)
    store_data(os.path.join(data_dir, 'smresidual.dat'), smresidual0)
    store_data(os.path.join(data_dir, 'ssmmpsf.dat'), ssmmpsf)
    store_data(os.path.join(data_dir, 'hsmmpsf.dat'), hsmmpsf)
    store_data(os.path.join(data_dir, 'ihsmmpsf.dat'), ihsmmpsf)
    store_data(os.path.join(data_dir, 'ldirty.dat'), ldirty)
    store_data(os.path.join(data_dir, 'psf.dat'), psf)



if __name__ == '__main__':
    np.random.seed(0)
    data_dir = './data/'

    test_cleaners(data_dir)
