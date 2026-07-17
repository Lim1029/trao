
# this code identifies signal window of a spectral cube, and save it as a hdu
# last update: 17 july 2026

# from useful_functions import baseline_cube
from spectral_cube import SpectralCube, BooleanArrayMask
import numpy as np
import astropy.units as u
from pybaselines.polynomial import imodpoly
import matplotlib.pyplot as plt
u.add_enabled_units(u.def_unit(['K (Tmb)'], represents=u.K))
u.add_enabled_units(u.def_unit(["K (Ta*)"], represents=u.K))
from codes.utils import cube_spectral_smooth, spectrum_smooth
from codes.utils.cube_mask import cube_mask
from codes.utils.spectrum_smooth import spectrum_smooth
from codes.utils.average_plotting import average_plotting
from pybaselines import Baseline
import sys 
from astropy.io import fits
from codes.utils.cube_smooth_tophat import TopHat_3DFilter
from scipy.ndimage import gaussian_filter1d, binary_dilation

# matplotlib configuration'
plt.style.use('./codes/astro.mplstyle')
plt.rcParams['figure.figsize'] = (1920/162, 1080/162)

############################ cube operation ############################
def cube_trim(cube, vmin, vmax):
    cube = cube.spectral_slab(vmin*u.km/u.s, vmax*u.km/u.s)
    return cube
    
def cube_smooth(cube, target_reso):
    cube_smoothed_spectral = cube_spectral_smooth(cube, target_reso, unit=u.km/u.s)
    cube = cube_smoothed_spectral
    return cube


def jcmt_window_test(cube, avg_mode='mean'):

    raw_data = cube.filled_data[:].value # Keep a pristine backup of the raw data
    data = raw_data.copy()

    smooth_kernel=[5,5,5]
    clips=[3,3,3]
    
    data = TopHat_3DFilter(data, smooth_kernel)
        
    nv, ny, nx = data.shape
    n_iters = len(clips)

    vel = cube.spectral_axis.to(u.km/u.s).value
    dv = np.abs(vel[1] - vel[0])  # Velocity spacing per channel in km/s

    # fixing the bin size to 5 km/s ( km/s to number of channels)
    bin_size = int(np.round(5.0 / dv))
    if bin_size < 1:
        bin_size = 1  

    #Calculate how many bins fit into the total channels
    nbin = nv // bin_size

    print('Number of bins: ', nbin)

    remainder = nv % bin_size 
    
    edges = [j * bin_size for j in range(nbin + 1)]
    if remainder != 0:
        edges[-1] = nv

    # Track the mask layout across each discrete clip iteration (4D array)
    window = np.ones((n_iters, nv, ny, nx), dtype=bool)
    
    for x in range(nx):
        for y in range(ny):
            spec = data[:, y, x]
            if np.all(np.isnan(spec)):
                continue
            
            # Create a mutable copy to aggressively wipe out peaks step-by-step
            spec_working = spec.copy()
            
            for i, clip in enumerate(clips):
                # Compute binned averages based on the working spectrum
                binned = [np.nanmedian(spec_working[edges[j]:edges[j+1]]) for j in range(nbin)]
                binned = np.array(binned, dtype=float)
                
                if np.all(np.isnan(binned)):
                    break

                if avg_mode == 'mean':
                    avg = np.nanmean(binned)
                elif avg_mode == 'median':
                    avg = np.nanmedian(binned)

                # 2. Calculate the absolute deviations from the median
                med_binned = np.nanmedian(binned)
                abs_deviation = np.abs(binned - med_binned)

                # 3. Get the median of those deviations
                mad = np.nanmedian(abs_deviation)

                # 4. Scale it to match standard deviation (if MAD is 0 due to flat data, fallback to std)
                if mad > 0:
                    std = 1.4826 * mad
                else:
                    std = np.nanstd(binned)

                binned_clipped = np.where(binned > avg + clip * std, np.nan, binned)
                mask_binned = ~np.isnan(binned_clipped)
                
                # --- MAP BACK TO FULL CHANNEL SCALE ---
                if remainder == 0: 
                    mask_full = np.repeat(mask_binned, bin_size)
                else: 
                    mask_full = np.repeat(mask_binned[:-1], bin_size)
                    mask_full = np.append(mask_full, [mask_binned[-1] for _ in range(bin_size + remainder)])
                mask_full = np.array(mask_full)
                
                # --- 10-CHANNEL WING BUFFER ---
                emission_full = ~mask_full
                
                if np.sum(emission_full) > 0: 
                    emission_expanded = binary_dilation(emission_full, structure=np.ones(11))
                    mask_full = ~emission_expanded
                
                # --- BANDPASS EDGE PROTECTION ---
                mask_full[:5] = True
                mask_full[-5:] = True
                
                # Clean out current features so they don't corrupt the next pass' noise calculations
                spec_working[~mask_full] = np.nan
                window[i, :, y, x] = mask_full
                
    # Keeps only channels marked True (noise) across ALL iterations
    master_window = np.logical_and.reduce(window, axis=0) 
                
    fig, axes = plt.subplots(3, 3, figsize=(12, 10), sharex=True)
    axes = axes.flatten()
    
    cube_max = np.nanmax(data, axis=0)
    max_sort = np.argsort(cube_max.flatten())
    nspec_total = len(max_sort)
    
    top3 = max_sort[-4:-1]
    mid3 = max_sort[nspec_total // 2 : nspec_total // 2 + 3]
    low3 = max_sort[:3]
    all_chosen = np.concatenate([top3, mid3, low3])
    all_chosen = np.unravel_index(all_chosen, cube_max.shape)
    
    mask_colors = ['red', 'magenta', 'cyan', 'orange', 'green']
    
    for i in range(9):
        iy, ix = all_chosen[0][i], all_chosen[1][i]
        ax = axes[i]
        
        spec_smoothed = data[:, iy, ix]
        spec_raw = raw_data[:, iy, ix] # Pull the backup raw spec
    
        # Plot raw data in light gray in the background, and smoothed data on top
        ax.step(vel, spec_raw, color='lightgray', linewidth=0.5, alpha=0.5, label='Raw Spec')
        ax.step(vel, spec_smoothed, color='black', linewidth=0.8, alpha=0.7, label='Smoothed Spec')
        
        for iter_idx in range(n_iters):
            window_iter = window[iter_idx, :, iy, ix]
            color = mask_colors[iter_idx % len(mask_colors)]
            mask_display = window_iter * np.nanmax(spec_smoothed) * 0.9
            ax.step(vel, mask_display, color=color, linewidth=1.2, 
                    alpha=0.8, label=fr"Iter {iter_idx + 1} ({clips[iter_idx]}$\sigma$, bins={nbin})") 
        ax.set_title(f"(x={ix}, y={iy})", fontsize=8)
        if i == 0:  
            ax.legend(fontsize=7, loc='upper right')
            
    plt.tight_layout()
    plt.show()
    
    return master_window, edges, nbin

if __name__ == "__main__":

    comment = 'During the signal windowing process, '
    cube_path = input("Input full/relative path to cube fits: ")
    cube = SpectralCube.read(cube_path).with_spectral_unit(u.km/u.s)
   
    print(f"Currently, the cube has a bandwidth of {cube.spectral_axis[0]:.2f} to {cube.spectral_axis[-1]:.2f} km/s")
    vlims = input("Input vmin and vmax in km/s separated by ',' to trim the cube (0 to skip trimming): ")
    if vlims != '0':
        vmin, vmax = list(map(float,vlims.split(',')))
        cube = cube_trim(cube, vmin, vmax)
        comment = comment + f"The cube is trimmed to [{vlims}] km/s. "

    print("masking....")
    window, edges, nbin = jcmt_window_test(cube)


    comment = comment + f"Window is automatically defined with nbin of {nbin} with clips [3,3,3]. "
   
    show_plot = input("Visualise baseline fitting? (y/n): ").strip().lower()

    if show_plot == 'y':
    
        snr_threshold = float(input("Enter SNR threshold for cube masking (0 to skip): "))

        average_plotting(
            cube,
            window=window,
            snr_threshold=snr_threshold
        )
    
    
    # saving results by creating a file with multiple hdu
    output_path = input("input the output path (enter n to skip saving): ")
    if output_path != 'n':
        primary_header = cube.header
        primary_header.add_comment(comment)
        hdu1 = fits.PrimaryHDU(data=cube.unmasked_data[:].value, header=cube.header)
        hdu2 = fits.ImageHDU(data=window, header=primary_header, name='WINDOW')
        hdul = fits.HDUList([hdu1, hdu2])
        hdul.writeto(output_path, overwrite=True)
        print(f'file saved to {output_path}')   
