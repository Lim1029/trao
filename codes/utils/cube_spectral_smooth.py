from astropy.convolution import Gaussian1DKernel
import astropy.units as u
import numpy as np
from spectral_cube import SpectralCube

def cube_spectral_smooth(cube, target_reso, unit=u.km/u.s, resample=True):
    fwhm_to_sigma = np.sqrt(8*np.log(2))
    spectral_axis = cube.spectral_axis.to(unit)
    current_reso = np.abs(np.diff(spectral_axis)[0].value)
    fwhm_gaussian = np.sqrt(target_reso**2 - current_reso**2)
    stddev_gaussian = fwhm_gaussian / fwhm_to_sigma 
    stddev_gaussian_pixel = stddev_gaussian/current_reso
    # we want the kernel in pixel units, so we force to km/s and take the value
    spectral_smoothing_kernel = Gaussian1DKernel(stddev=stddev_gaussian_pixel)
    cube_spectralsmoothed = cube.spectral_smooth(spectral_smoothing_kernel)
    # account for reverse velocity axis (e.g., IRAM)
    if resample:
        if np.diff(spectral_axis)[0].value < 0:
            new_spectral_axis = np.arange(spectral_axis[0].value, spectral_axis[-1].value, -target_reso)
        else:
            new_spectral_axis = np.arange(spectral_axis[0].value, spectral_axis[-1].value, target_reso)
        cube_spectralsmoothed_spectralresample = cube_spectralsmoothed.spectral_interpolate(new_spectral_axis*unit,suppress_smooth_warning=True)
        return cube_spectralsmoothed_spectralresample
    else:
        return cube_spectralsmoothed
 
if __name__ == "__main__":
    cube_path = input("input cube path: ")
    cube = SpectralCube.read(cube_path)
    target_reso = float(input("input target reso (km/s): "))
    cube_smoothed = cube_spectral_smooth(cube, target_reso)
    print("done smoothing")
    output_path = input("input output path to save the cube: ")
    cube_smoothed.write(output_path, overwrite=True)
