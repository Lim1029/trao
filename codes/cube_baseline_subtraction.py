# this code explores different methods of performing baseline subtraction of a cube
# last update: 8 May 2026

# from useful_functions import baseline_cube
from spectral_cube import SpectralCube, BooleanArrayMask
import numpy as np
import astropy.units as u
from pybaselines.polynomial import imodpoly
import matplotlib.pyplot as plt
u.add_enabled_units(u.def_unit(['K (Tmb)'], represents=u.K))
u.add_enabled_units(u.def_unit(["K (Ta*)"], represents=u.K))
# cube = SpectralCube.read(args.input_cube_path)
# D:/NaritNextcloud/NARIT_RA/carta/TRAO/fits_files/w43
from utils.cube_spectral_smooth import cube_spectral_smooth
from utils.average_plotting import average_plotting
from utils.cube_mask import cube_mask
from utils.spectrum_smooth import spectrum_smooth
from pybaselines import Baseline
import sys 
from astropy.io import fits
from cube_signal_window import jcmt_window
import matplotlib.colors as mcolors
from scipy.interpolate import interp1d
from scipy.interpolate import make_splrep, make_lsq_spline
import math

############################ cube operation ############################
def cube_trim(cube, vmin, vmax):
    cube = cube.spectral_slab(vmin*u.km/u.s, vmax*u.km/u.s)
    return cube
    
def cube_smooth(cube, target_reso):
    cube_smoothed_spectral = cube_spectral_smooth(cube, target_reso, unit=u.km/u.s)
    return cube_smoothed_spectral

############################ fitting method ############################
def baseline_pixel_poly(cube, window, poly_order):
    x = cube.spectral_axis.value
    data = cube.unmasked_data[:].value # (nspec, ny, nx)
    mask_arr = window
    nspec, ny, nx = data.shape
    corrected_cube = np.empty_like(data)
    baseline_cube = np.empty_like(data)
    for j in range(ny):
        for i in range(nx):
            y = data[:,j,i]
            mask = mask_arr[:,j,i]
            baseline = imodpoly(y, x, poly_order=poly_order, weights=mask,
                mask_initial_peaks=False, num_std=3, use_original=True)[0]
            corrected = y - baseline
            corrected_cube[:,j,i] = corrected
            baseline_cube[:,j,i] = baseline
    return corrected_cube, baseline_cube
# this code shall iteratively perform baseline fitting with different poly order, and stop at the lowest AIC

# adopted from the class code by Slawa Kabanovic
def baseline_iterative_poly(cube, window, max_order, smooth_vel=0.3):
    x = cube.spectral_axis.value
    data = cube.unmasked_data[:].value # (nspec, ny, nx)
    mask_arr = window
    nspec, ny, nx = data.shape
    corrected_cube = np.empty_like(data)
    baseline_cube = np.empty_like(data)
    polyorder_cube = np.empty_like(data[0])
    
    if smooth_vel != 0:
        cube_smoothed = cube_spectral_smooth(cube, smooth_vel)
        data_smoothed = cube_smoothed.unmasked_data[:].value
        x_smoothed = cube_smoothed.spectral_axis.value
    else:
        data_smoothed = data
        x_smoothed = x
    
    for j in range(ny):
        for i in range(nx):
            y = data_smoothed[:,j,i]
            y_ori = data[:,j,i]
            mask = mask_arr[:,j,i]
            #downsample the mask
            target_length = len(y)
            current_indices = np.linspace(0, len(mask) - 1, len(mask))
            target_indices = np.linspace(0, len(mask) - 1, target_length)
            f = interp1d(current_indices, mask, kind='nearest')
            mask = f(target_indices)
            
            if np.any(np.isnan(y)) or np.any(np.isnan(mask)):
                continue
            aic_old = 9999
            is_continue = 1
            poly_order = 0
            selected_polyorder = 0
            # breakpoint()
            while is_continue and (poly_order <= max_order):
                
                baseline = imodpoly(y, x_smoothed, poly_order=poly_order, weights=mask,
                                    mask_initial_peaks=False, num_std=1, use_original=True)[0]
                # breakpoint()
                # calculate aic 
                baselined = y - baseline #aka residual
                n = sum(mask)
                rss = np.sum((baselined * mask)**2)
                k = poly_order + 1
                aic = n * np.log(rss/n) + 2*k
                # print(n,rss,k,aic)
                
                # check if the fit is 'doing too much', by checking the derivation, 
                # only check if this is not the first round
                if poly_order != 0:                       
                    sigma = np.std(baselined)
                    signal_indices = np.where(mask == 0)[0]
                    # print(signal_indices)
                    derivative = np.full(len(y), np.nan)
                    vel_res = np.diff(x)[0]
                    
                    for k in signal_indices:
                        if k == 0:
                            derivative[k] = (baseline[k+1] - baseline[k]) / vel_res
                        elif k == len(y) - 1:
                            derivative[k] = (baseline[k] - baseline[k-1]) / vel_res
                        else:
                            derivative[k] = (baseline[k+1] - baseline[k-1]) / (2.0 * vel_res)

                    find_extrema = False
                    
                    # check only interior signal indices
                    for k in signal_indices:
                        if k == 0:
                            continue
                        if mask[k-1] != 0:
                            continue

                        # sign change in derivative
                        if derivative[k] * derivative[k-1] < 0:

                            noise_limit = 2.0 * sigma
                            poly_diff = baseline[k]-baseline_old[k]

                            if abs(poly_diff) > noise_limit:
                                find_extrema = True
                                # print('found extrema!')
                                selected_polyorder = poly_order - 1
                                is_continue = 0
                                break
                                
                if aic >= aic_old:
                    # use the previous poly_order instead
                    selected_polyorder = poly_order - 1
                    is_continue = 0
                    
                poly_order += 1
                aic_old = aic
                baseline_old = baseline
                
            # print(f'best fit is order {selected_polyorder}')
                
            baseline,params = imodpoly(y, x_smoothed, poly_order=selected_polyorder, weights=mask,
                                mask_initial_peaks=False, num_std=3, use_original=True, return_coef=True)
            coef = params['coef']
            # apply the baseline to the cube of original shape 
            baseline_full = np.polynomial.Polynomial(coef)(x)
            corrected = y_ori - baseline_full
            corrected_cube[:,j,i] = corrected
            baseline_cube[:,j,i] = baseline_full
            polyorder_cube[j,i] = poly_order

    # plotting of the selected polyorder
    cmap = plt.get_cmap('viridis', max_order+1)
    bounds = np.arange(-0.5, max_order+1.5, 1)
    norm = mcolors.BoundaryNorm(bounds, cmap.N)
    fig, ax = plt.subplots(figsize=(8, 6))
    im = ax.imshow(polyorder_cube, cmap=cmap, norm=norm)
    cbar = plt.colorbar(im, ticks=np.arange(0, max_order+1))
    cbar.set_label('Polynomial Order')
    plt.show()
            
    return corrected_cube, baseline_cube
    
# reference: https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.make_splrep.html#scipy.interpolate.make_splrep    
def baseline_spline(cube, window,knot_start= 100, knot_spacing = 100, k=3):
    x = cube.spectral_axis.value
    data = cube.unmasked_data[:].value # (nspec, ny, nx)
    mask_arr = window
    nspec, ny, nx = data.shape
    corrected_cube = np.empty_like(data)
    baseline_cube = np.empty_like(data)
    for j in range(ny):
        for i in range(nx):
            y = data[:,j,i]
            y = np.nan_to_num(y)
            mask = mask_arr[:,j,i]
            mask = mask.ravel()
            # make knots (t) every 100 channels, and at the boundary, but skip the emission window
            t = x[np.where(mask==1)][knot_start:-1:knot_spacing]
            t = np.r_[(x[0],)*4,t,(x[-1],)*4]       

            # breakpoint()
            try:
                spl = make_lsq_spline(x, y, t, k, w=mask)
            except:
                breakpoint()
            baseline = spl(x)
            
            corrected = y - baseline
            corrected_cube[:,j,i] = corrected
            baseline_cube[:,j,i] = baseline
    return corrected_cube, baseline_cube    
############################ window defining ############################
# def jcmt_window(cube, nbin):
    # window, edges = jcmt_window(cube, nbin=nbin)
    # return window
    
def fix_window(cube, wlist):
    x = cube.spectral_axis.value
    total_window = np.full(len(x), True)
    for idx in range(0,len(wlist),2):
        wleft, wright = wlist[idx], wlist[idx+1]
        window = ((x < wleft) | (x > wright))
        total_window &= window
    window = np.broadcast_to(total_window[:, None, None], cube.shape).astype(int)
    return window
    
def no_window(cube):
    x = cube.spectral_axis.value
    total_window = np.full(len(x), True)
    window = np.broadcast_to(total_window[:, None, None], cube.shape).astype(int)
    return window  
    
    
if __name__ == "__main__":
    comment = 'During the baseline subtraction process, '
    cube_path = input("Input full/relative path to cube fits: ")
    cube = SpectralCube.read(cube_path).with_spectral_unit(u.km/u.s)
    print(f"Currently, the cube has a bandwidth of {cube.spectral_axis[0]:.2f} to {cube.spectral_axis[-1]:.2f} km/s")
    vlims = input("Input vmin and vmax in km/s separated by ',' to trim the cube (0 to skip trimming): ")
    if vlims != '0':
        vmin, vmax = list(map(float,vlims.split(',')))
        cube = cube_trim(cube, vmin, vmax)
        comment = comment + f"The cube is trimmed to [{vlims}] km/s. "
    target_reso = float(input("Input target resolution (km/s) to perform smoothing prior to baseline fitting (0 to skip): "))
    if target_reso != 0:
        cube = cube_smooth(cube, target_reso=target_reso)
        comment = comment + f"The cube is smoothed spectrally to {target_reso} km/s. "
    window_method = int(input("How to define spectral window? (1) no window (2) user-defined (3) automatic: "))
    match window_method:
        case 1:
            window = no_window(cube)
            comment = comment + "No window is defined. "
        case 2:
            vlims = input("Input window pair list separated by ',' (e.g., 0,15,25,40) [km/s]: ")
            vlims = list(map(float,vlims.split(',')))
            window = fix_window(cube, vlims)
            comment = comment + f'Fix spectral window of [{vlims}] is defined. '
        case 3:
            nbin = int(input(f"Input nbin, to divide {len(cube.spectral_axis)} channels: "))
            clips = input("Input clip (in unit of sigma) separated by , (e.g., 2,2.5,3): ")
            clips = list(map(float,clips.split(',')))
            smooth_kernel = input("tophat smooth cube prior to windowing? (enter kernel size v,y,x or 0 to skip): ")
            if smooth_kernel != '0':
                smooth_kernel = list(map(int,smooth_kernel.split(',')))
            else:
                smooth_kernel = None
            window, edges = jcmt_window(cube, nbin=nbin, clips=clips, smooth_kernel=smooth_kernel)
            comment = comment + f"Window is automatically defined with nbin of {nbin}, clipping with {clips} and tophat smoothed with kernel size {smooth_kernel} prior to windowing."
    fitting_method = int(input("Input method (1) fixed poly (2) iterative poly (3) spline: "))
    match fitting_method:
        case 1:
            poly_order = int(input("Input poly order: "))
            corrected, baseline = baseline_pixel_poly(cube, window, poly_order)
            comment = comment + f"The baseline is fitted with polynomial function with order {poly_order}. "
        case 2:
            max_order = int(input("Input maximum poly order to iterate: "))
            corrected, baseline = baseline_iterative_poly(cube, window, max_order)
            comment = comment + f"The baseline is fitted with polynomial function with automatically selected order, until order {max_order}. "
        case 3:
            knot_start = int(input("Input knot starting channel:"))
            knot_spacing = int(input("Input knot spacing in channels:"))
            corrected, baseline = baseline_spline(cube, window, knot_start = knot_start, knot_spacing = knot_spacing)
            comment = comment + f"The baseline is fitted with cube spline function."
        
        case _:
            print('Nothing is done')
            sys.exit(0)
            
            
    # visualising the baseline fitting
    show_plot = input(
        "Visualise baseline fitting? (y/n): "
    ).strip().lower()

    if show_plot == 'y':
        
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
    hdu4 = fits.ImageHDU(data=window, header=cube.header, name='MASK')
    hdul = fits.HDUList([hdu1, hdu2, hdu3, hdu4])
    hdul.writeto(output_path, overwrite=True)
