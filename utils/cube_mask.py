from codes.cube_rms_estimate import rms_negative
import numpy as np

# this mask is to make any PIXEL with snr < threshold become nan
def cube_mask(cube, mask_threshold):
    spectra = cube.unmasked_data[:].value
    spectra_max = np.nanmax(spectra, axis=0)
    spectra_rms,_ = rms_negative(cube)
    spectra_snr = spectra_max/spectra_rms
    mask = spectra_snr>mask_threshold

    expanded_mask = mask[np.newaxis, :, :]
    masked_spectra = np.where(expanded_mask, spectra, np.nan)
    return masked_spectra
    
# this mask is to make any CHANNEL with value < threshold of that PIXEL become nan
def cube_mask_channel(cube, mask_threshold):
    spectra = cube.unmasked_data[:].value
    spectra_rms,_ = rms_negative(cube)
    # Force the 2D RMS map to match the 3D cube shape
    rms_3d = np.broadcast_to(spectra_rms, spectra.shape)
    mask = spectra > (mask_threshold * rms_3d)
    masked_spectra = np.where(mask, spectra, np.nan)
    return masked_spectra