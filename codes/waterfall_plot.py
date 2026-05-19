# this code plots a waterfall plot similar to Fig. 6 in Higgins+21
# the x axis is velocity, the y axis is the flatten pixel number (we will have ~3600), then do a imshow.

from spectral_cube import SpectralCube
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from codes.utils import cube_spectral_smooth
plt.style.use('./codes/astro.mplstyle')

def waterfall_plot(cube, ax):
    # flatten the x y axis of the cube
    spectra = cube.unmasked_data[:].value
    x_axis = cube.spectral_axis.value
    nv, ny, nx = spectra.shape
    spectra_flattened = spectra.reshape(nv, ny*nx).T
    vmin = np.nanpercentile(spectra_flattened, 20)
    vmax = np.nanpercentile(spectra_flattened, 99)
    im = ax.imshow(spectra_flattened, origin='lower', vmin=vmin, vmax=vmax, aspect='auto')
    
    # ax.set_xticks(x_axis)
    return im, ax
    
    
if __name__ == '__main__':
    cube_path = input("input full/relative path to cube fits: ")
    cube = SpectralCube.read(cube_path).with_spectral_unit(u.km/u.s)
    target_reso = float(input("Input target resolution (km/s) to perform smoothing prior to plotting (0 to skip): "))
    if target_reso != 0:
        cube = cube_spectral_smooth(cube, target_reso, unit=u.km/u.s)
    fig = plt.figure(figsize=(4,2))
    ax = fig.add_subplot(111)
    im,ax = waterfall_plot(cube, ax)
    plt.colorbar(im, ax=ax)
    plt.show()