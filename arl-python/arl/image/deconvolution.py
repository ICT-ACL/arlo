""" Image deconvolution functions

The standard deconvolution algorithms are provided:

    hogbom: Hogbom CLEAN See: Hogbom CLEAN A&A Suppl, 15, 417, (1974)
    
    msclean: MultiScale CLEAN See: Cornwell, T.J., Multiscale CLEAN (IEEE Journal of Selected Topics in Sig Proc,
    2008 vol. 2 pp. 793-801)

    mfsmsclean: MultiScale Multi-Frequency See: U. Rau and T. J. Cornwell, “A multi-scale multi-frequency
    deconvolution algorithm for synthesis imaging in radio interferometry,” A&A 532, A71 (2011).

For example to make dirty image and PSF, deconvolve, and then restore::

    model = create_image_from_visibility(vt, cellsize=0.001, npixel=256)
    dirty, sumwt = invert_2d(vt, model)
    psf, sumwt = invert_2d(vt, model, dopsf=True)

    comp, residual = deconvolve_cube(dirty, psf, niter=1000, threshold=0.001, fracthresh=0.01, window='quarter',
                                 gain=0.7, algorithm='msclean', scales=[0, 3, 10, 30])

    restored = restore_cube(comp, psf, residual)

"""

import numpy
import logging

from astropy.convolution import Gaussian2DKernel, convolve
from photutils import fit_2dgaussian

from arl.data.data_models import Image
from arl.data.parameters import get_parameter
from arl.image.operations import create_image_from_array, copy_image, \
    calculate_image_frequency_moments, calculate_image_from_frequency_moments

from arl.image.cleaners import hogbom, msclean, msmfsclean

log = logging.getLogger(__name__)

def deconvolve_cube(dirty: Image, psf: Image, **kwargs) -> (Image, Image):
    """ Clean using a variety of algorithms
    
    Functions that clean a dirty image using a point spread function. The algorithms available are:
    
    hogbom: Hogbom CLEAN See: Hogbom CLEAN A&A Suppl, 15, 417, (1974)
    
    msclean: MultiScale CLEAN See: Cornwell, T.J., Multiscale CLEAN (IEEE Journal of Selected Topics in Sig Proc,
    2008 vol. 2 pp. 793-801)

    mfsmsclean: MultiScale Multi-Frequency See: U. Rau and T. J. Cornwell, “A multi-scale multi-frequency
    deconvolution algorithm for synthesis imaging in radio interferometry,” A&A 532, A71 (2011).
    
    For example::
    
        comp, residual = deconvolve_cube(dirty, psf, niter=1000, gain=0.7, algorithm='msclean',
                                         scales=[0, 3, 10, 30], threshold=0.01)
    
    :param dirty: Image dirty image
    :param psf: Image Point Spread Function
    :param window: Window image (Bool) - clean where True
    :param algorithm: Cleaning algorithm: 'msclean'|'hogbom'|'mfsmsclean'
    :param gain: loop gain (float) 0.7
    :param threshold: Clean threshold (0.0)
    :param fracthres: Fractional threshold (0.01)
    :param scales: Scales (in pixels) for multiscale ([0, 3, 10, 30])
    :param nmoments: Number of frequency moments (default 3)
    :param findpeak: Method of finding peak in mfsclean: 'Algorithm1'|'ASKAPSoft'|'CASA'|'ARL', Default is ARL.
    :return: componentimage, residual
    
    """
    assert type(dirty) == Image, "Type is %s" % (type(dirty))
    assert type(psf) == Image, "Type is %s" % (type(psf))

    opt = get_parameter(kwargs, 'opt', False)
    if not opt:
        log.debug('Using original algorithm')
    else:
        log.debug('Using optimized algorithm')
    
    window = get_parameter(kwargs, 'window', None)
    if window == 'quarter':
        qx = dirty.shape[3] // 4
        qy = dirty.shape[2] // 4
        window = numpy.zeros_like(dirty.data)
        window[..., (qy + 1):3 * qy, (qx + 1):3 * qx] = 1.0
        log.info('deconvolve_cube: Cleaning inner quarter of each sky plane')
    else:
        window = None
    
    psf_support = get_parameter(kwargs, 'psf_support', None)
    if isinstance(psf_support, int):
        if (psf_support < psf.shape[2] // 2) and ((psf_support < psf.shape[3] // 2)):
            centre = [psf.shape[2] // 2, psf.shape[3] // 2]
            psf.data = psf.data[..., (centre[0] - psf_support):(centre[0] + psf_support),
                       (centre[1] - psf_support):(centre[1] + psf_support)]
            log.info('deconvolve_cube: PSF support = +/- %d pixels' % (psf_support))
    
    algorithm = get_parameter(kwargs, 'algorithm', 'msclean')
    
    if algorithm == 'msclean':
        log.info("deconvolve_cube: Multi-scale clean of each polarisation and channel separately")
        gain = get_parameter(kwargs, 'gain', 0.7)
        assert 0.0 < gain < 2.0, "Loop gain must be between 0 and 2"
        thresh = get_parameter(kwargs, 'threshold', 0.0)
        assert thresh >= 0.0
        niter = get_parameter(kwargs, 'niter', 100)
        assert niter > 0
        scales = get_parameter(kwargs, 'scales', [0, 3, 10, 30])
        fracthresh = get_parameter(kwargs, 'fractional_threshold', 0.01)
        assert 0.0 < fracthresh < 1.0
    
        comp_array = numpy.zeros_like(dirty.data)
        residual_array = numpy.zeros_like(dirty.data)
        for channel in range(dirty.data.shape[0]):
            for pol in range(dirty.data.shape[1]):
                if psf.data[channel, pol, :, :].max():
                    log.info("deconvolve_cube: Processing pol %d, channel %d" % (pol, channel))
                    if window is None:
                        comp_array[channel, pol, :, :], residual_array[channel, pol, :, :] = \
                            msclean(dirty.data[channel, pol, :, :], psf.data[channel, pol, :, :],
                                    None, gain, thresh, niter, scales, fracthresh)
                    else:
                        comp_array[channel, pol, :, :], residual_array[channel, pol, :, :] = \
                            msclean(dirty.data[channel, pol, :, :], psf.data[channel, pol, :, :],
                                    window[channel, pol, :, :], gain, thresh, niter, scales, fracthresh)
                else:
                    log.info("deconvolve_cube: Skipping pol %d, channel %d" % (pol, channel))
                    
        comp_image = create_image_from_array(comp_array, dirty.wcs)
        residual_image = create_image_from_array(residual_array, dirty.wcs)

    elif algorithm == 'msmfsclean' or algorithm == 'mfsmsclean':
        findpeak = get_parameter(kwargs, "findpeak", 'ARL')
        
        log.info("deconvolve_cube: Multi-scale multi-frequency clean of each polarisation separately")
        nmoments = get_parameter(kwargs, "nmoments", 3)
        assert nmoments > 0, "Number of frequency moments must be greater than zero"
        dirty_taylor = calculate_image_frequency_moments(dirty, nmoments=nmoments)
        psf_taylor = calculate_image_frequency_moments(psf, nmoments=2 * nmoments)
    
        gain = get_parameter(kwargs, 'gain', 0.7)
        assert 0.0 < gain < 2.0, "Loop gain must be between 0 and 2"
        thresh = get_parameter(kwargs, 'threshold', 0.0)
        assert thresh >= 0.0
        niter = get_parameter(kwargs, 'niter', 100)
        assert niter > 0
        scales = get_parameter(kwargs, 'scales', [0, 3, 10, 30])
        fracthresh = get_parameter(kwargs, 'fractional_threshold', 0.1)
        assert 0.0 < fracthresh < 1.0

        comp_array = numpy.zeros(dirty_taylor.data.shape)
        residual_array = numpy.zeros(dirty_taylor.data.shape)
        for pol in range(dirty_taylor.data.shape[1]):
            if psf_taylor.data[0, pol, :, :].max():
                log.info("deconvolve_cube: Processing pol %d" % (pol))
                if window is None:
                    comp_array[:, pol, :, :], residual_array[:, pol, :, :] = \
                        msmfsclean(dirty_taylor.data[:, pol, :, :], psf_taylor.data[:, pol, :, :],
                                None, gain, thresh, niter, scales, fracthresh, findpeak, opt)
                else:
                    comp_array[:, pol, :, :], residual_array[:, pol, :, :] = \
                        msmfsclean(dirty_taylor.data[:, pol, :, :], psf_taylor.data[:, pol, :, :],
                                window[:, pol, :, :], gain, thresh, niter, scales, fracthresh, findpeak, opt)
            else:
                log.info("deconvolve_cube: Skipping pol %d" % (pol))

        comp_image = create_image_from_array(comp_array, dirty_taylor.wcs)
        residual_image = create_image_from_array(residual_array, dirty_taylor.wcs)
        
        return_moments = get_parameter(kwargs, "return_moments", False)
        if not return_moments:
            log.info("Deconvolve_cube: calculating spectral cubes")
            comp_image = calculate_image_from_frequency_moments(dirty, comp_image)
            residual_image = calculate_image_from_frequency_moments(dirty, residual_image)
        else:
            log.info("Deconvolve_cube: constructed moment cubes")

    elif algorithm == 'hogbom':
        log.info("deconvolve_cube: Hogbom clean of each polarisation and channel separately")
        gain = get_parameter(kwargs, 'gain', 0.7)
        assert 0.0 < gain < 2.0, "Loop gain must be between 0 and 2"
        thresh = get_parameter(kwargs, 'threshold', 0.0)
        assert thresh >= 0.0
        niter = get_parameter(kwargs, 'niter', 100)
        assert niter > 0
        fracthresh = get_parameter(kwargs, 'fractional_threshold', 0.1)
        assert 0.0 <= fracthresh < 1.0
        
        comp_array = numpy.zeros(dirty.data.shape)
        residual_array = numpy.zeros(dirty.data.shape)
        for channel in range(dirty.data.shape[0]):
            for pol in range(dirty.data.shape[1]):
                if psf.data[channel, pol, :, :].max():
                    log.info("deconvolve_cube: Processing pol %d, channel %d" % (pol, channel))
                    if window is None:
                        comp_array[channel, pol, :, :], residual_array[channel, pol, :, :] = \
                            hogbom(dirty.data[channel, pol, :, :], psf.data[channel, pol, :, :],
                                   None, gain, thresh, niter, fracthresh)
                    else:
                        comp_array[channel, pol, :, :], residual_array[channel, pol, :, :] = \
                            hogbom(dirty.data[channel, pol, :, :], psf.data[channel, pol, :, :],
                                   window[channel, pol, :, :], gain, thresh, niter, fracthresh)
                else:
                    log.info("deconvolve_cube: Skipping pol %d, channel %d" % (pol, channel))
        
        comp_image = create_image_from_array(comp_array, dirty.wcs)
        residual_image = create_image_from_array(residual_array, dirty.wcs)
    else:
        raise ValueError('deconvolve_cube: Unknown algorithm %s' % algorithm)
    
    return comp_image, residual_image


def restore_cube(model: Image, psf: Image, residual=None, **kwargs) -> Image:
    """ Restore the model image to the residuals

    :params psf: Input PSF
    :return: restored image

    """
    assert type(model) == Image, "Type is %s" % (type(model))
    assert type(psf) == Image, "Type is %s" % (type(psf))
    assert residual is None or type(residual) == Image, "Type is %s" % (type(residual))

    restored = copy_image(model)
    
    npixel = psf.data.shape[3]
    sl = slice(npixel // 2 - 7, npixel // 2 + 8)
    
    # isotropic at the moment!
    try:
        fit = fit_2dgaussian(psf.data[0, 0, sl, sl])
        if fit.x_stddev <= 0.0 or fit.y_stddev <= 0.0:
            log.debug('restore_cube: error in fitting to psf, using 1 pixel stddev')
            size = 1.0
        else:
            size = max(fit.x_stddev, fit.y_stddev)
            log.debug('restore_cube: psfwidth = %s' % (size))
    except:
        log.debug('restore_cube: warning in fit to psf, using 1 pixel stddev')
        size = 1.0
    
    # By convention, we normalise the peak not the integral so this is the volume of the Gaussian
    norm = 2.0 * numpy.pi * size ** 2
    gk = Gaussian2DKernel(size)
    for chan in range(model.shape[0]):
        for pol in range(model.shape[1]):
            restored.data[chan, pol, :, :] = norm * convolve(model.data[chan, pol, :, :], gk, normalize_kernel=False)
    if residual is not None:
        restored.data += residual.data
    return restored
