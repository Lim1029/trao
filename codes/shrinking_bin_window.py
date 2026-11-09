# this code identifies signal window of a spectral cube, and saves it as an HDU

import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d, binary_dilation
from pybaselines import Baseline
from pybaselines.polynomial import imodpoly

import astropy.units as u
from astropy.io import fits
from spectral_cube import SpectralCube, BooleanArrayMask

# Custom module imports
u.add_enabled_units(u.def_unit(['K (Tmb)'], represents=u.K))
u.add_enabled_units(u.def_unit(["K (Ta*)"], represents=u.K))

from codes.utils import cube_spectral_smooth, spectrum_smooth
from codes.utils.cube_mask import cube_mask
from codes.utils.average_plotting import average_plotting
from codes.utils.cube_smooth_tophat import TopHat_3DFilter

# Matplotlib configuration
plt.style.use('./codes/astro.mplstyle')
plt.rcParams['figure.figsize'] = (1920 / 162, 1080 / 162)


# trimming the ends
def cube_trim(cube, vmin, vmax):
    return cube.spectral_slab(vmin * u.km / u.s, vmax * u.km / u.s)

#smoothing 
def cube_smooth(cube, target_reso):
    return cube_spectral_smooth(cube, target_reso, unit=u.km / u.s)

# windowing 
def jcmt_window_test(cube, avg_mode='mean'):

    raw_data = cube.filled_data[:].value  # Keeping the raw data 
    data = raw_data.copy()

    smooth_kernel = [5, 5, 5]
    clips = [3, 3, 3]

    data = TopHat_3DFilter(data, smooth_kernel)

    nv, ny, nx = data.shape
    n_iters = len(clips)

    vel = cube.spectral_axis.to(u.km / u.s).value
    dv = np.abs(vel[1] - vel[0])  # Velocity spacing per channel

    # Base bin size to ~5 km/s
    base_bin_size = int(np.round(5.0 / dv))
    if base_bin_size < 1:
        base_bin_size = 1

    base_nbin = nv // base_bin_size
    print('Initial number of base bins:', base_nbin)

    # Track 4D mask layout across iterations
    window = np.ones((n_iters, nv, ny, nx), dtype=bool)

    edges = None # placeholder for return reference
    for x in range(nx):
        for y in range(ny):
            spec = data[:, y, x]
            
            if np.all(np.isnan(spec)):
                continue

            spec_working = spec.copy()

            for i, clip in enumerate(clips):
            
                # Dynamically calculate current bin parameters per iteration (shrinking/growing multi-scale)
                current_nbin = base_nbin * (2**i)
                if current_nbin > nv:
                    current_nbin = nv

                # Split spectrum into exact contiguous chunks to guarantee total length matches nv precisely
                chunks = np.array_split(spec_working, current_nbin)
                binned = np.array([np.nanmean(chunk) for chunk in chunks], dtype=float)

                if np.all(np.isnan(binned)):
                    break

                if avg_mode == 'mean':
                    avg = np.nanmean(binned)
                elif avg_mode == 'median':
                    avg = np.nanmedian(binned)

                # Absolute deviations & MAD
                med_binned = np.nanmedian(binned)
                abs_deviation = np.abs(binned - med_binned)
                mad = np.nanmedian(abs_deviation)

                std = 1.4826 * mad if mad > 0 else np.nanstd(binned)

                binned_clipped = np.where(binned > avg + clip * std, np.nan, binned)
                mask_binned = ~np.isnan(binned_clipped)

                # Map back to channel scale precisely matching each chunk's length
                mask_full = np.concatenate([np.repeat(mb, len(chunk)) for mb, chunk in zip(mask_binned, chunks)])
                mask_full = np.array(mask_full, dtype=bool)

                # Signal wing dilation
                emission_full = ~mask_full
                if np.sum(emission_full) > 0:
                    emission_expanded = binary_dilation(emission_full, structure=np.ones(11))
                    mask_full = ~emission_expanded

                # Bandpass edge protection
                mask_full[:5] = True
                mask_full[-5:] = True

                # Exclude found features for subsequent passes
                spec_working[~mask_full] = np.nan
                window[i, :, y, x] = mask_full

    # Combine noise regions across iterations
    master_window = np.logical_and.reduce(window, axis=0)

    # Visualization grid
    fig, axes = plt.subplots(3, 3, figsize=(12, 10), sharex=True)
    axes = axes.flatten()

    cube_max = np.nanmax(data, axis=0)
    max_sort = np.argsort(cube_max.flatten())
    nspec_total = len(max_sort)

    top3 = max_sort[-4:-1]
    mid3 = max_sort[nspec_total // 2: nspec_total // 2 + 3]
    low3 = max_sort[:3]
    all_chosen = np.unravel_index(np.concatenate([top3, mid3, low3]), cube_max.shape)

    mask_colors = ['red', 'magenta', 'cyan', 'orange', 'green']

    for i in range(9):
        iy, ix = all_chosen[0][i], all_chosen[1][i]
        ax = axes[i]

        spec_smoothed = data[:, iy, ix]
        spec_raw = raw_data[:, iy, ix]

        ax.step(vel, spec_raw, color='lightgray', linewidth=0.5, alpha=0.5, label='Raw Spec')
        ax.step(vel, spec_smoothed, color='black', linewidth=0.8, alpha=0.7, label='Smoothed Spec')

        for iter_idx in range(n_iters):
            window_iter = window[iter_idx, :, iy, ix]
            color = mask_colors[iter_idx % len(mask_colors)]
            mask_display = window_iter * np.nanmax(spec_smoothed) * 0.9
            ax.step(vel, mask_display, color=color, linewidth=1.2,
                    alpha=0.8, label=rf"Iter {iter_idx + 1} ({clips[iter_idx]} sigma)")
        ax.set_title(f"(x={ix}, y={iy})", fontsize=8)
        if i == 0:
            ax.legend(fontsize=7, loc='upper right')

    plt.tight_layout()
    plt.show()

    return master_window, base_nbin, clips, smooth_kernel


if __name__ == "__main__":

    comment = 'During the signal windowing process, '
    cube_path = input("Input full/relative path to cube fits: ")
    cube = SpectralCube.read(cube_path).with_spectral_unit(u.km / u.s)

    print(f"Currently, the cube has a bandwidth of {cube.spectral_axis[0]:.2f} to {cube.spectral_axis[-1]:.2f} km/s")
    vlims = input("Input vmin and vmax in km/s separated by ',' to trim the cube (0 to skip trimming): ")
    if vlims != '0':
        vmin, vmax = list(map(float, vlims.split(',')))
        cube = cube_trim(cube, vmin, vmax)
        comment += f"The cube is trimmed to [{vlims}] km/s. "

    snr_threshold = float(input('Mask the cube based on SNR? (enter 0 to skip masking): '))
    if snr_threshold != 0:
        from codes.cube_rms_estimate import rms_negative

        spectra_max = np.max(cube, axis=0).value
        spectra_rms, _ = rms_negative(cube)
        spectra_snr = spectra_max / spectra_rms

        mask = spectra_snr > snr_threshold
        expanded_mask = mask[np.newaxis, :, :]
        cube = cube.with_mask(expanded_mask)

    window, nbin, clips, smooth_kernel = jcmt_window_test(cube)

    comment += f"Window is automatically defined with multi-scale iterative binning and clips [3,3,3]. "

    show_plot = input("Visualise baseline fitting? (y/n): ").strip().lower()
    if show_plot == 'y':
        snr_threshold = float(input("Enter SNR threshold for cube masking (0 to skip): "))
        average_plotting(
            cube,
            window=window,
            snr_threshold=snr_threshold
        )

    # Save output to FITS file
    output_path = input("Input the output path (enter n to skip saving): ").strip()
    if output_path != 'n':
        primary_header = cube.header.copy()
        primary_header.add_comment(comment)
        hdu1 = fits.PrimaryHDU(data=cube.unmasked_data[:].value, header=cube.header)
        hdu2 = fits.ImageHDU(data=window, header=primary_header, name='WINDOW')
        hdul = fits.HDUList([hdu1, hdu2])
        hdul.writeto(output_path, overwrite=True)
        print(f'File saved to {output_path}')
