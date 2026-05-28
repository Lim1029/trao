# this codes plot the average spectrum of raw cube and show the masking window additionally the average corrected spectrum, and baseline can be plotted
# masking using signal to noise ratio threshold is also an option.
# last updated: 8 may 2026


import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u

from codes.utils.cube_mask import cube_mask
from codes.utils.spectrum_smooth import spectrum_smooth

def average_plotting(cube, corrected=None, window=None, baseline=None, snr_threshold=0):


    # cube: raw spectra 
     
    vel = cube.spectral_axis.to(u.km/u.s).value
    
    # averaging and masking(optional) the corrected spectra, baseline and raw spectra.
    
    # no masking 
    if snr_threshold == 0:
    
        original_avg = cube.mean(axis=(1,2)).value
        
        if corrected is not None:
            corrected_avg  = np.nanmean(corrected, axis=(1,2))
        
        if baseline is not None:
            baseline_avg  = np.nanmean(baseline, axis=(1,2))
        
       
    # with masking 
    
    else:
        original_masked, expanded_mask = cube_mask(cube, snr_threshold) 
        original_avg = np.nanmean(original_masked, axis=(1,2))
        
        if corrected is not None:
            corrected_masked = np.where(expanded_mask, corrected, np.nan)
            corrected_avg = np.nanmean(corrected_masked, axis=(1,2))
            
        if baseline is not None:
            baseline_masked = np.where(expanded_mask, baseline, np.nan)
            baseline_avg = np.nanmean(baseline_masked, axis=(1,2))
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.step(vel, original_avg, linewidth=0.5, alpha=0.2, color='black', label='raw')
  
    if baseline is not None:
        ax.plot(vel, baseline_avg, color= 'black', linestyle='--', label='baseline')
    
    if corrected is not None:
        averaged_smoothed = spectrum_smooth(corrected_avg, vel, 0.3)
        ax.step(averaged_smoothed.spectral_axis, averaged_smoothed.flux, label='baselined_smoothed')
            
    
    original_averaged_smoothed = spectrum_smooth(original_avg, vel, 0.3)
    ax.step(original_averaged_smoothed.spectral_axis, original_averaged_smoothed.flux, alpha=0.7, color='orange', label='raw_smoothed')
    
    if window is not None:
        window_averaged = np.nanmean(window, axis=(1,2))
        
        
        for i in range(len(vel)-1):
            ax.axvspan(vel[i], vel[i+1], alpha=(1-window_averaged[i]), color='yellow')
            
    
    ax.set_ylabel('K (Ta*)')
    ax.set_xlabel('LSR Velocity (km/s)')
    ax.axhline(0, color='red', linestyle='dashed')
    #ax.legend()

    #if snr_threshold == 0:
     #   ax.set_title("Average spectra, no masking")
    #else:
     #   ax.set_title(f"Average spectra, mask SNR > {snr_threshold}")
    plt.show()
    
    # saving option for plot
    
    save_plot = input(
        "Save this plot? (y/n): "
    ).strip().lower()

    if save_plot == 'y':

        save_path = input(
            "Enter output path to save plot: "
        )

        fig.savefig(
            save_path,
            dpi=300,
            bbox_inches='tight'
        )

        print(f"Plot saved to {save_path}")
