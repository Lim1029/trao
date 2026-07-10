from specutils import Spectrum1D
import numpy as np
from astropy.convolution import Gaussian1DKernel, convolve
from specutils.manipulation import gaussian_smooth, SplineInterpolatedResampler
import astropy.units as u

# unitless for all parameters
def spectrum_smooth(spectrum, spectral_axis, target_resolution, resample=True):
    if spectral_axis[1]-spectral_axis[0] < 0: # to account for reverse spectral axis issue, like IRAM cube
        input_spec = Spectrum1D(spectral_axis=spectral_axis[::-1]*(u.km/u.s), flux=spectrum[::-1]*u.K)
        new_spectral_axis = np.arange(spectral_axis[-1], spectral_axis[0], target_resolution)
    else:
        input_spec = Spectrum1D(spectral_axis=spectral_axis*(u.km/u.s), flux=spectrum*u.K)
        new_spectral_axis = np.arange(spectral_axis[0], spectral_axis[-1], target_resolution)
    
    current_resolution = abs(spectral_axis[1]-spectral_axis[0])
    fwhm_gaussian = (target_resolution**2 - current_resolution**2)**0.5
    sigma_gaussian = fwhm_gaussian / (2*np.sqrt(2*np.log(2)))
    sigma_gaussian_pixel = sigma_gaussian / current_resolution
    # print(f"Number of pixels of the kernel: {sigma_gaussian_pixel}")    
    spec1_gsmooth = gaussian_smooth(input_spec, stddev=sigma_gaussian_pixel)
    
    kernel = Gaussian1DKernel(stddev=sigma_gaussian_pixel)
    smoothed_spectrum = convolve(spectrum, kernel)
<<<<<<< Updated upstream
    
    resampler = SplineInterpolatedResampler()
    # resampler = FluxConservingResampler()
    # resampler = LinearInterpolatedResampler()
    new_spectrum = resampler(spec1_gsmooth, new_spectral_axis*(u.km/u.s))
=======
    if resample:
        resampler = SplineInterpolatedResampler()
        # resampler = FluxConservingResampler()
        # resampler = LinearInterpolatedResampler()
        new_spectrum = resampler(spec1_gsmooth, new_spectral_axis*(u.km/u.s))
    else:
        new_spectrum = smoothed_spectrum
>>>>>>> Stashed changes
    return new_spectrum
