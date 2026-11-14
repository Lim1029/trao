# This codes is written to model a common baseline for the entire cube using the 
# those spatial pixels which are assigned as completely emission free (no emission at all along spectra axis.)
# NOTE- the selected pixels are subject to masking.
#
# the median value is used here to model the baseline. For the baseline modelling polynomial function is used
# 

import sys 
import numpy as np
import astropy.units as u
from astropy.io import fits
import matplotlib.pyplot as plt
from spectral_cube import SpectralCube
from pybaselines.polynomial import imodpoly
from codes.shrinking_bin_window import jcmt_window_test
from codes.utils.average_plotting import average_plotting

u.add_enabled_units(u.def_unit(['K (Tmb)'], represents=u.K))
u.add_enabled_units(u.def_unit(["K (Ta*)"], represents=u.K))

def common_baseline_white_pixels(cube, window, poly_order):

    x = cube.spectral_axis.value
    data = cube.unmasked_data[:].value  # Shape: (nspec, ny, nx)
    nspec, ny, nx = data.shape
    
    # Find spatial pixels that are line-free across ALL channels
    # window shape is (nspec, ny, nx). Checking along axis=0 gives a (ny, nx) boolean mask.

    pure_bkg_mask = np.all(window > 0, axis=0)
    num_clean_pixels = np.sum(pure_bkg_mask)
    print(f"Found {num_clean_pixels} fully clean(white pixels) spatial pixels out of {ny * nx} total pixels.")
    
    # Safety check for data starvation
    if num_clean_pixels == 0:
        raise ValueError(
            "Error: No spatial pixels are completely line-free across all channels! "
        )
    
    clean_spectra = data[:, pure_bkg_mask]  # Shape: (nspec, num_clean_pixels)
    
    med_spec = np.nanmedian(clean_spectra, axis=1)

    #breakpoint()
    
    weights = np.ones(nspec, dtype=int)
    baseline_1d, params = imodpoly(
        med_spec, x, poly_order=poly_order, weights=weights,
        mask_initial_peaks=False, num_std=3, use_original=True, return_coef=True
    )
    coef = params['coef']
    
    fig, axes = plt.subplots(2, 1, figsize=(10, 10))
    
    # Plot to check the location of the white pixels.
    im = axes[0].imshow(pure_bkg_mask, origin='lower', cmap='gray')
    axes[0].set_title(f"Pure Bkg Mask ({num_clean_pixels} pixels)")
    axes[0].set_xlabel("X Pixel")
    axes[0].set_ylabel("Y Pixel")
    
    # Plot of the median spectrum and the polynomial fit
    axes[1].plot(x, med_spec, label='Median spectra', color='gray', alpha=0.7)
    axes[1].plot(x, baseline_1d, label=f'Poly Fit (order={poly_order})', color='red', lw=2)
    axes[1].axhline(y=0, ls='--', color='black', lw=0.7)
    axes[1].set_title("Median spectra & Baseline Fit")
    axes[1].set_xlabel("Spectral Axis (km/s)")
    axes[1].set_ylabel("Intensity")
    axes[1].legend()
    
    fig.subplots_adjust(hspace=0.3)
    
    plt.show()

    # Broadcast 1D baseline back to 3D cube shape and subtract

    baseline_cube = np.broadcast_to(baseline_1d[:, None, None], data.shape).copy()
    corrected_cube = data - baseline_cube
    
    return corrected_cube, baseline_cube



if __name__ == "__main__":

    comment = 'During the baseline subtraction process, '
    cube_path = input("Input full/relative path to cube fits: ")
    cube = SpectralCube.read(cube_path).with_spectral_unit(u.km/u.s)

    # Drop the first and last 20 spectral channels
    cube = cube[20:-20, :, :]
    #comment = comment + "First and last 20 channels only contains NaN values. "

    print(f"Currently, the cube has a bandwidth of {cube.spectral_axis[0]:.2f} to {cube.spectral_axis[-1]:.2f} km/s")

    print("Windowing with JCMT method with shrinking bin (initial bin size 5 km/s)...")
    window, nbin, clips, smooth_kernel = jcmt_window_test(cube)
    comment = comment + f"Window is automatically defined with nbin of {nbin}, clipping with {clips} and tophat smoothed with kernel size {smooth_kernel} prior to windowing."
            
    print("Baselining with a common baseline estimated using pixels which are completely emission free across all spectral channels...")
    poly_order = int(input("Input poly order:"))
    corrected, baseline = common_baseline_white_pixels(cube, window, poly_order)
    comment = comment + f"The baseline is fitted with polynomial function with order {poly_order} and common baseline modeled from averaged spectra. "
   

    # visualising the baseline fitting
    print("Visualise baseline fitting... ")

    snr_threshold = float(input("Enter SNR threshold for cube masking (0 to skip): "))

    average_plotting(
            cube,
            corrected=corrected,
            window=window,
            baseline=baseline,
            snr_threshold=snr_threshold
    )       
    # saving results by creating a file with multiple hdu
    output_path = input("input the output path: ")
    primary_header = cube.header
    primary_header.add_comment(comment)
    hdu1 = fits.PrimaryHDU(data=corrected, header=primary_header)
    hdu2 = fits.ImageHDU(data=baseline, header=cube.header, name='BASELINE')
    hdu3 = fits.ImageHDU(data=cube.unmasked_data[:].value, header=cube.header, name='UN-BASELINED')
    hdu4 = fits.ImageHDU(data=window.astype(int), header=cube.header, name='MASK')
    hdul = fits.HDUList([hdu1, hdu2, hdu3, hdu4])
    hdul.writeto(output_path, overwrite=True)
