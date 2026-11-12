# this code identifies signal window of a spectral cube, and save it as a hdu
# last update: 24 March 2026

# from useful_functions import baseline_cube
from spectral_cube import SpectralCube, BooleanArrayMask
import numpy as np
import astropy.units as u
from pybaselines.polynomial import imodpoly
import matplotlib.pyplot as plt
u.add_enabled_units(u.def_unit(['K (Tmb)'], represents=u.K))
u.add_enabled_units(u.def_unit(["K (Ta*)"], represents=u.K))
from utils import cube_spectral_smooth, spectrum_smooth
from utils.cube_mask import cube_mask
from utils.spectrum_smooth import spectrum_smooth
from pybaselines import Baseline
import sys 
from astropy.io import fits
from utils.cube_smooth_tophat import TopHat_3DFilter
from scipy.ndimage import gaussian_filter1d, binary_dilation

# matplotlib configuration'
plt.style.use('./astro.mplstyle')
plt.rcParams['figure.figsize'] = (1920/162, 1080/162)

############################ cube operation ############################
def cube_trim(cube, vmin, vmax):
    cube = cube.spectral_slab(vmin*u.km/u.s, vmax*u.km/u.s)
    return cube
    
def cube_smooth(cube, target_reso):
    cube_smoothed_spectral = cube_spectral_smooth(cube, target_reso, unit=u.km/u.s)
    cube = cube_smoothed_spectral
    return cube

############################ window defining ############################
# this follows the MFITTREND procedure to identify spectral window dynamically

def jcmt_window_spectrum(v,spec,nbin=20, clips=[2,2.5,3], avg_mode='mean'):
    nv = len(v)
    bin_size = nv // nbin
    remainder = nv % nbin 
    edges = [i * bin_size for i in range(nbin+1)]
    if remainder != 0:
        edges[-1] = nv
    binned = [np.mean(spec[edges[i]:edges[i+1]]) for i in range(nbin)]
    binned = np.array(binned, dtype=float)
    # iteratively calculate mean and standard deviation of binned means and mask out outliers 
    for clip in clips:
        if avg_mode == 'mean':
            avg = np.nanmean(binned)
        elif avg_mode == 'median':
            avg = np.nanmedian(binned)
        std = np.nanstd(binned)
        binned = np.where(binned > avg+clip*std, np.nan, binned)
                
    mask_binned = ~np.isnan(binned)
    # expand the masked binned back to its original shape
    if remainder == 0: 
        mask_full = np.repeat(mask_binned, bin_size)
    else: # for non dividible nbin, the last bin would have size of remainder
        mask_full = np.repeat(mask_binned[:-1], bin_size)
        mask_full = np.append(mask_full, [mask_binned[-1] for i in range(bin_size+remainder)])
    mask_full = np.array(mask_full)
    return mask_full


def browse_profiles(data, vel, window=None, nrows=3, ncols=3, descending=True, window_on_twin_axis=True):

    nperp = nrows * ncols

    # ---- sort every pixel by its peak value -------------------------------
    cube_max = np.nanmax(data, axis=0)              # (ny, nx)
    flat_max = cube_max.ravel()
    order    = np.argsort(flat_max)                 # ascending
    order    = order[np.isfinite(flat_max[order])]  # drop NaN pixels
    if descending:
        order = order[::-1]                         # strongest first
    iy_all, ix_all = np.unravel_index(order, cube_max.shape)

    nspec = order.size
    npages = max(1, int(np.ceil(nspec / nperp)))
    state = {'page': 0}

    # ---- figure -----------------------------------------------------------
    fig, axes = plt.subplots(nrows, ncols, figsize=(12, 8), squeeze=False)
    axes = axes.flatten()
    fig.subplots_adjust(hspace=0.5, wspace=0.3)
    twins = [None] * nperp

    def draw(page):
        page = int(np.clip(page, 0, npages - 1))
        state['page'] = page
        start = page * nperp
        stop  = min(start + nperp, nspec)

        for k, ax in enumerate(axes):
            j = start + k
            ax.clear()
            if twins[k] is not None:
                twins[k].remove(); twins[k] = None
            if j >= stop:                      # nothing left on this page
                ax.set_visible(False)
                continue
            ax.set_visible(True)

            iy, ix = iy_all[j], ix_all[j]
            ax.step(vel, data[:, iy, ix], where='mid', lw=0.9, color='C0')

            if window is not None:
                w = np.asarray(window[:, iy, ix])
                if window_on_twin_axis:
                    ax2 = ax.twinx()
                    ax2.step(vel, w, where='mid', lw=0.8,
                             color='C1', alpha=0.7)
                    ax2.tick_params(labelsize=6, colors='C1')
                    twins[k] = ax2
                else:
                    ax.step(vel, w, where='mid', lw=0.8,
                            color='C1', alpha=0.7)

            ax.set_title(f"x={ix}, y={iy}  |  peak={flat_max[order[j]]:.3g}",
                         fontsize=8)
            ax.tick_params(labelsize=7)

        fig.suptitle(
            f"{'strongest → weakest' if descending else 'weakest → strongest'}   "
            f"page {page+1}/{npages}   "
            f"(profiles {start+1}–{stop} of {nspec})\n"
            "[n]ext  [p]rev  [home]first  [end]last  [s]ave  [q]uit",
            fontsize=9)
        fig.canvas.draw_idle()

    def on_key(event):
        key = (event.key or '').lower()
        if   key in ('n', 'right'): draw(state['page'] + 1)
        elif key in ('p', 'left'):  draw(state['page'] - 1)
        elif key == 'home':         draw(0)
        elif key == 'end':          draw(npages - 1)
        elif key == 's':
            fig.savefig(f"profiles_page_{state['page']+1:03d}.png", dpi=150)
        elif key in ('q', 'escape'):
            plt.close(fig)

    fig.canvas.mpl_connect('key_press_event', on_key)
    draw(0)
    plt.show()
    return fig


def jcmt_window(cube, nbin=30, clips=[2,2.5,3], plot_progress=None, avg_mode='mean', smooth_kernel=None, bin_expand=1):
    data = cube.filled_data[:].value # (nspec, ny, nx)
    if smooth_kernel != None:   # for weak low SNR emission, we may want to smooth first
        data = TopHat_3DFilter(data, smooth_kernel)
    window = np.full(data.shape, np.nan)
    vel = cube.spectral_axis.to(u.km/u.s).value
    nv, ny, nx = data.shape # define the binning edges
    bin_size = nv // nbin
    remainder = nv % nbin # we will get remainder if the choice of nbin is not divisible by nv
    edges = [i * bin_size for i in range(nbin+1)] # define the [) inclusive exclusive edges of bins
    if remainder != 0: # if have remainder, we put them in the last bin
        edges[-1] = nv
    
    for x in range(nx):
        for y in range(ny):
            spec = data[:,y,x]
            if np.all(np.isnan(spec)):
                continue
            # print(spec)
            # compute the mean value of each bin
            binned = [np.mean(spec[edges[i]:edges[i+1]]) for i in range(nbin)]
            binned = np.array(binned, dtype=float)
            
            # iteratively calculate mean and standard deviation of binned means and mask out outliers 
            for clip in clips:
                if avg_mode == 'mean':
                    avg = np.nanmean(binned)
                elif avg_mode == 'median':
                    avg = np.nanmedian(binned)
                std = np.nanstd(binned)
                binned = np.where(binned > avg+clip*std, np.nan, binned)
                
                # we can plot and see
                if plot_progress:
                    plt.step(vel, spec, linewidth=0.5, alpha=0.5)
                    for edge in edges[:-1]:
                        plt.axvline(vel[edge], linestyle='dashed', linewidth=0.5, alpha=0.5, color='red')
                    if remainder == 0:
                        binned_expand = np.repeat(binned, bin_size)
                    else:
                        binned_expand = np.repeat(binned, bin_size)
                        binned_expand = np.append(binned_expand, [binned[-1] for i in range(remainder)])

                    plt.step(vel, binned_expand, linewidth=1)
                    plt.axhline(avg, linewidth=0.5, color='black')
                    plt.axhline(avg+std, linewidth=0.5, color='orange')
                    plt.axhline(avg+clip*std, linewidth=1, color='red')
                    plt.show()
                
            mask_binned = ~np.isnan(binned)
            mask_binned_old = mask_binned.copy()

            # protect the emission wings, by assigning true to the neighour
            emission_binned = ~mask_binned
            emission_binned = binary_dilation(emission_binned, iterations=bin_expand)
            mask_binned = ~emission_binned
            
            # protect the boundary to exclude from masking (important during spline baseline fitting)
            mask_binned[0] = True
            mask_binned[-1] = True
            # expand the masked binned back to its original shape
            if remainder == 0: 
                mask_full = np.repeat(mask_binned, bin_size)
            else: # for non dividible nbin, the last bin would have size of remainder
                mask_full = np.repeat(mask_binned[:-1], bin_size)
                mask_full = np.append(mask_full, [mask_binned[-1] for i in range(bin_size+remainder)])
            mask_full = np.array(mask_full)
            window[:,y,x] = mask_full
        
    browse_profiles(data, vel, window)
    
    return window, edges

def window_visualise(cube, window,baseline, original, snr_threshold=0):

    vel = cube.spectral_axis.value  
    
    # averaging the corrected spectra, baseline and raw spectra.
    # no masking 
    if snr_threshold == 0:
        corrected_avg = cube.mean(axis=(1,2))
        baseline_avg  = np.nanmean(baseline, axis=(1,2))
        original_avg  = np.nanmean(original, axis=(1,2))
       
    # with masking 
    else:
        corrected_masked, expanded_mask = cube_mask(cube, snr_threshold) 
        baseline_masked = np.where(expanded_mask, baseline, np.nan)
        original_masked = np.where(expanded_mask, original, np.nan)

        corrected_avg = np.nanmean(corrected_masked, axis=(1,2))
        baseline_avg = np.nanmean(baseline_masked, axis=(1,2))
        original_avg = np.nanmean(original_masked, axis=(1,2))
    
    window_averaged = np.nanmean(window, axis=(1,2))
    
    # print(window_averaged)
    # is_window = np.where(window_averaged>0)
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.step(vel, original_avg, linewidth=0.5, alpha=0.1, color='black')
    ax.step(vel, corrected_avg, linewidth=0.5, alpha=0.1, color='blue')
    ax.plot(vel, baseline_avg, linewidth = 0.5, color= 'black', linestyle='--')
    averaged_smoothed = spectrum_smooth(corrected_avg, vel, 0.3)
    original_averaged_smoothed = spectrum_smooth(original_avg, vel, 0.3)
    # breakpoint()
    ax.step(averaged_smoothed.spectral_axis, averaged_smoothed.flux)
    ax.step(original_averaged_smoothed.spectral_axis, original_averaged_smoothed.flux, alpha=0.5)
    ax.set_ylabel('K (T$_{A}$)')
    ax.set_xlabel('LSR Velocity (km/s)')
    for i in range(len(vel)-1):
        ax.axvspan(vel[i], vel[i+1], alpha=1-window_averaged[i], color='yellow')
    ax.axhline(0, color='red', linestyle='dashed')
    # ax.step(vel, window_averaged)
    # for edge in edges[:-1]:
        # ax.axvline(vel[edge], linestyle='dashed', linewidth=0.5, alpha=0.5, color='red')
        # also write the window fraction
        # ax.text(x=vel[edge],y=ax.get_ylim()[1],s=f"{1-window_averaged[edge]:.2f}",ha='left', va='top')
    # if snr_threshold == 0:
        # ax.set_title("Average spectra, no masking")
    # else:
        # ax.set_title(f"Average spectra, mask SNR > {snr_threshold}")
    plt.show()

        
if __name__ == "__main__":
    choice = input('plot window fraction (1) or perform dynamic window (2)?: ')
    if choice == '1':
        cube_path = input("Input full/relative path to cube fits: ")
        cube = SpectralCube.read(cube_path).with_spectral_unit(u.km/u.s)
        snr_threshold = float(input("Enter SNR threshold for cube masking (0 to skip): "))
        hdul = fits.open(cube_path)
        window = hdul['mask'].data
        original = hdul['UN-BASELINED'].data
        baseline = hdul['BASELINE'].data
        window_visualise(cube, window, baseline, original, snr_threshold)
    elif choice == '2':
        comment = 'During the signal windowing process, '
        cube_path = input("Input full/relative path to cube fits: ")
        cube = SpectralCube.read(cube_path).with_spectral_unit(u.km/u.s)
        print(f"Currently, the cube has a bandwidth of {cube.spectral_axis[0]:.2f} to {cube.spectral_axis[-1]:.2f} km/s")
        vlims = input("Input vmin and vmax in km/s separated by ',' to trim the cube (0 to skip trimming): ")
        if vlims != '0':
            vmin, vmax = list(map(float,vlims.split(',')))
            cube = cube_trim(cube, vmin, vmax)
            comment = comment + f"The cube is trimmed to [{vlims}] km/s. "
        
        snr_threshold = float(input('mask the cube based on SNR? (enter 0 to skip masking): '))
        if snr_threshold != 0:
            from codes.cube_rms_estimate import rms_negative
            spectra_max = np.max(cube, axis=0).value
            spectra_rms,_ = rms_negative(cube)
            spectra_snr = spectra_max/spectra_rms
            # breakpoint()
            mask = spectra_snr > snr_threshold
            expanded_mask = mask[np.newaxis, :, :]
            cube = cube.with_mask(expanded_mask)
            
        print('performing dynamic singal windowing...')
        nbin = int(input(f"Input nbin, to divide {len(cube.spectral_axis)} channels: "))
        clips = input("Input clip (in unit of sigma) separated by , (e.g., 2,2.5,3): ")
        clips = list(map(float,clips.split(',')))
        to_show = input("show the windowing process (y/n)?: ")
        if to_show == 'y':
            to_show = True
        else:
            to_show = False
        smooth_kernel = input("tophat smooth cube prior to windowing? (enter kernel size v,y,x or 0 to skip): ")
        if smooth_kernel != '0':
            smooth_kernel = list(map(int,smooth_kernel.split(',')))
        else:
            smooth_kernel = None
        bin_expand = int(input("Number of neighboring bin each side to mask as emission to protect wings (default 1): "))
        window, edges = jcmt_window(cube, nbin=nbin, clips=clips, smooth_kernel=smooth_kernel, bin_expand= bin_expand, plot_progress=to_show)
        comment = comment + f"Window is automatically defined with nbin of {nbin} with clips {clips}. "
                    
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
            
    else:
        print('do nothing')
        sys.exit(0)
