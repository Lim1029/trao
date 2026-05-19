# this is an extension to the typical cube_baseline_subtraction,
# to be used in individual OTF of TRAO W40 C18O and CN 
# dynamic signal windowing is first determined on the averaged cube
# and then applied to individual OTFs to perform baseline fitting
# individual baselined-OTF is then averaged again to obtain an averaged cube.
# the reason not to identify signal window on individual OTF is because of the very weak SNR (tested and failed)

# for C18O
# otfs_folder_path = 'carta/TRAO/jihye/W40/C18O/fits_files/'
# otf_filenames = [f'W40_C18O_0{i}.fits' for i in range(3,10)]
# otfs_output_folder = 'carta/TRAO/jihye/W40/C18O/fits_files/baselined/'
# cube_path = 'carta/TRAO/fits_files/w40/c18o.fits'
# output_cube_path = 'carta/TRAO/jihye/W40/C18O/fits_files/baselined/spline/W40_C18O_03-09.fits'
# vmin, vmax = -50,50
# nbin=20
# clips=[2,2.5,3]

# for CN
otfs_folder_path = 'carta/TRAO/jihye/W40/CN/fits_files/'
otf_filenames = [f'W40_CN_0{i}.fits' for i in range(3,10)]
otfs_output_folder = 'carta/TRAO/jihye/W40/CN/fits_files/baselined/'
cube_path = 'carta/TRAO/fits_files/w40/cn.fits'
output_cube_path = 'carta/TRAO/jihye/W40/CN/fits_files/baselined/spline/W40_CN_03-09.fits'
vmin, vmax = -50,50
nbin=15
clips=[2,2.5,3]

from spectral_cube import SpectralCube
import astropy.units as u
from astropy.io import fits
from codes.cube_rms_estimate import rms_negative

import numpy as np

cube = SpectralCube.read(cube_path).spectral_slab(vmin*u.km/u.s,vmax*u.km/u.s)

from codes.cube_signal_window import jcmt_window

window, edges = jcmt_window(cube, nbin=nbin, clips=clips)

from codes.cube_baseline_subtraction import baseline_iterative_poly, baseline_spline

baselined_otfs = []
array_rms = []

import matplotlib.pyplot as plt

for filename in otf_filenames:
    print(f'performing otf {filename}')
    otf = SpectralCube.read(otfs_folder_path+filename).spectral_slab(vmin*u.km/u.s,vmax*u.km/u.s)
    # corrected, baseline = baseline_iterative_poly(otf, window, 5)
    corrected, baseline = baseline_spline(otf, window)
    # breakpoint()
    baselined_otfs.append(corrected)
    rms_map, mean_rms = rms_negative(corrected)
    # plt.imshow(rms_map, origin='lower')
    # plt.colorbar()
    # plt.show()
    array_rms.append(rms_map)

baselined_otfs = np.array(baselined_otfs)
array_rms = np.array(array_rms)
# calculate the weights of each pixel in each scan 
array_weight = 1/array_rms**2

# calculate the simple or weighted average of all scans
weights_expanded = array_weight[:, np.newaxis, :, :] 
weighted_mean_cube = np.sum(baselined_otfs * weights_expanded, axis=0) / np.sum(weights_expanded, axis=0)

# calculate and display the noise rms of the mean cube
weighted_mean_cube_rms, mean_rms = rms_negative(weighted_mean_cube)
print(mean_rms)

primary_header = cube.header
hdu1 = fits.PrimaryHDU(data=weighted_mean_cube, header=primary_header)
hdu2 = fits.ImageHDU(data=cube.unmasked_data[:].value, header=cube.header, name='UN-BASELINED')
hdu3 = fits.ImageHDU(data=window, header=cube.header, name='MASK')
hdul = fits.HDUList([hdu1,hdu2,hdu3])
hdul.writeto(output_cube_path, overwrite=True)

plt.imshow(weighted_mean_cube_rms, origin='lower')
plt.colorbar()
plt.show()