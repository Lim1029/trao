# this is an extension to the typical cube_baseline_subtraction,
# to be used in individual OTF of TRAO W40 C18O and CN 
# STEP 0:  First otf maps are combined using weighted averaging.
# STEP 1:  dynamic signal windowing is first determined on the averaged cube
# and then applied to individual OTFs to perform baseline fitting
# STEP 2:  individual baselined-OTF is then averaged again to obtain an averaged cube.
# STEP 3:  with this averaged cube again dynamic signal window is determined.
# STEP 4:  Again Step 2 is repeated to get a averaged cube.
# the reason not to identify signal window on individual OTF is because of the very weak SNR (tested and failed)

# for C18O
# otfs_folder_path = '/home/vinay/fits_files/c18o/'
# otf_filenames = [f'W40_C18O_0{i}.fits' for i in range(3,10)]
# otfs_output_folder = '/home/vinay/fits_files/c18o/baselined1/'
# cube_path = '/home/vinay/fits_files/c18o/w40/c18o.fits'
# output_cube_path = '/home/vinay/fits_files/c18o/baselined1/spline/W40_C18O_03-09.fits'
# vmin, vmax = -50,50
# nbin=20
# clips=[2,2.5,3]

# for CN
otfs_folder_path = '/home/vinay/fits_files/cn/'
otf_filenames = [f'W40_CN_0{i}.fits' for i in range(3,10)]
otfs_output_folder = '/home/vinay/fits_files/cn/baselined1/'  
cube_path = '/home/vinay/fits_files/cn/w40/cn.fits'
output_cube_path = '/home/vinay/fits_files/cn/baselined1/spline/W40_CN_03-09.fits'
vmin, vmax = -50,50
nbin=15
clips=[1.7,2,2.5]
smooth_kernel=[5,3,3]

from spectral_cube import SpectralCube
from astropy.wcs import WCS
import astropy.units as u
from astropy.io import fits
from codes.cube_rms_estimate import rms_negative
from codes.utils.average_plotting import average_plotting
from codes.cube_signal_window import jcmt_window
from codes.cube_baseline_subtraction import baseline_iterative_poly, baseline_spline
import numpy as np
import matplotlib.pyplot as plt

# reading the already combined cube for initial masking.
#cube = SpectralCube.read(cube_path).with_spectral_unit(u.km/u.s)
#cube = cube.spectral_slab(vmin*u.km/u.s,vmax*u.km/u.s)

reference_cube = (SpectralCube.read(otfs_folder_path + otf_filenames[0]).with_spectral_unit(u.km/u.s).spectral_slab(vmin*u.km/u.s,vmax*u.km/u.s))

primary_header = reference_cube.header
wcs = WCS(primary_header)

raw_otfs = []
array_rms0 = []

for filename in otf_filenames:

    otf = (SpectralCube.read(otfs_folder_path + filename).with_spectral_unit(u.km/u.s).spectral_slab(vmin*u.km/u.s,vmax*u.km/u.s))
    
    otf = otf.spectral_interpolate(reference_cube.spectral_axis)

    data = otf.unmasked_data[:].value
    
    print(filename)
    print(otf.shape)
    print(otf.spectral_axis[0])
    print(otf.spectral_axis[-1])
    print()

    raw_otfs.append(data)

    rms_map, mean_rms = rms_negative(data)

    array_rms0.append(rms_map)


raw_otfs = np.array(raw_otfs)
array_rms0 = np.array(array_rms0)

weights0 = 1 /(array_rms0)**2

weights_expanded0 = (weights0[:, np.newaxis, :, :])

combined_raw = (np.sum(raw_otfs * weights_expanded0,axis=0)/np.sum(weights_expanded0, axis=0))

combined_cube0 = SpectralCube(data=combined_raw * reference_cube.unit,wcs=wcs)

window, edges = jcmt_window(combined_cube0, nbin=nbin, clips=clips, smooth_kernel=smooth_kernel)


# First baselining 

baselined_otfs = []
array_rms = []

for filename in otf_filenames:
    print(f'performing otf {filename}')
    otf = SpectralCube.read(otfs_folder_path+filename).with_spectral_unit(u.km/u.s).spectral_slab(vmin*u.km/u.s,vmax*u.km/u.s)
    otf = otf.spectral_interpolate(reference_cube.spectral_axis)
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
primary_header = reference_cube.header
wcs =WCS(primary_header)

combined_cube = SpectralCube( data=weighted_mean_cube*reference_cube.unit, wcs=wcs)


# window masking using the combined cube with baselined otf maps

window2, edges2 = jcmt_window(combined_cube, nbin=nbin, clips=clips, smooth_kernel= smooth_kernel)


# second baselining

baselined_otfs2 = []
array_rms2 = []
baseline_otfs2 = []

for filename in otf_filenames:
    print(f'performing otf {filename}')
    otf = SpectralCube.read(otfs_folder_path+filename).with_spectral_unit(u.km/u.s).spectral_slab(vmin*u.km/u.s,vmax*u.km/u.s)
    otf = otf.spectral_interpolate(reference_cube.spectral_axis)
    # corrected, baseline = baseline_iterative_poly(otf, window, 5)
    corrected2, baseline2 = baseline_spline(otf, window2)
    # breakpoint()
    baselined_otfs2.append(corrected2)
    
    baseline_otfs2.append(baseline2)
    rms_map2, mean_rms2 = rms_negative(corrected2)
    # plt.imshow(rms_map, origin='lower')
    # plt.colorbar()
    # plt.show()
    array_rms2.append(rms_map2)
    
baselined_otfs2 = np.array(baselined_otfs2)
array_rms2 = np.array(array_rms2)
baseline_otfs2 = np.array(baseline_otfs2)

array_weight2 = 1 / array_rms2**2

weights_expanded2 = array_weight2[:, np.newaxis, :, :]

weighted_mean_cube2 = (
    np.sum(
        baselined_otfs2 * weights_expanded2,
        axis=0
    )
    /
    np.sum(weights_expanded2, axis=0)
)


weighted_baseline2 = (
    np.sum(
        baseline_otfs2 * weights_expanded2,
        axis=0
    )
    /
    np.sum(weights_expanded2, axis=0)
)

# calculate and display the noise rms of the mean cube
weighted_mean_cube_rms, mean_rms = rms_negative(weighted_mean_cube2)
print(mean_rms)


hdu1 = fits.PrimaryHDU(data=weighted_mean_cube2, header=primary_header)
hdu2 = fits.ImageHDU(data=combined_cube0.unmasked_data[:].value, header=reference_cube.header, name='UN-BASELINED')
hdu3 = fits.ImageHDU(data=window2, header=reference_cube.header, name='MASK')
hdul = fits.HDUList([hdu1,hdu2,hdu3])
hdul.writeto(output_cube_path, overwrite=True)

show_plot = input(
        "Visualise baseline fitting? (y/n): "
    ).strip().lower()

if show_plot == 'y':
        
        snr_threshold = float(input("Enter SNR threshold for cube masking (0 to skip): "))

        average_plotting(
            combined_cube0,
            corrected=weighted_mean_cube2,
            window=window2,
            baseline=weighted_baseline2,
            snr_threshold=snr_threshold
        )   

plt.imshow(weighted_mean_cube_rms, origin='lower')
plt.colorbar()
plt.show()
