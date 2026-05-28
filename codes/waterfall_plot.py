# this code plots a waterfall plot similar to Fig. 6 in Higgins+21
# the x axis is velocity, the y axis is the flatten pixel number (we will have ~3600), then do a imshow.

from spectral_cube import SpectralCube
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from codes.utils import cube_spectral_smooth
plt.style.use('./codes/astro.mplstyle')

def cube_trim(cube, vmin, vmax):
    cube = cube.spectral_slab(vmin*u.km/u.s, vmax*u.km/u.s)
    return cube

def waterfall_plot(cube, ax):
    # flatten the x y axis of the cube
    spectra = cube.unmasked_data[:].value
    x_axis = cube.spectral_axis.value
    nv, ny, nx = spectra.shape
    spectra_flattened = spectra.reshape(nv, ny*nx).T
    vmin = np.nanpercentile(spectra_flattened, 10)
    vmax = np.nanpercentile(spectra_flattened, 99)
    im = ax.imshow(spectra_flattened, origin='lower', vmin=vmin, vmax=vmax, aspect='auto', extent=[x_axis[0], x_axis[-1], 0, ny*nx] )
    
    # ax.set_xticks(x_axis)
    return im, ax
    
    
if __name__ == '__main__':
    cube_path = input("input full/relative path to cube fits: ")
    cube = SpectralCube.read(cube_path).with_spectral_unit(u.km/u.s)
    print(f"Currently, the cube has a bandwidth of {cube.spectral_axis[0]:.2f} to {cube.spectral_axis[-1]:.2f} km/s")
    vlims = input("Input vmin and vmax in km/s separated by ',' to trim the cube (0 to skip trimming): ")
    if vlims != '0':
        vmin, vmax = list(map(float,vlims.split(',')))
        cube = cube_trim(cube, vmin, vmax)
    target_reso = float(input("Input target resolution (km/s) to perform smoothing prior to plotting (0 to skip): "))
    if target_reso != 0:
        cube = cube_spectral_smooth(cube, target_reso, unit=u.km/u.s)
    fig = plt.figure(figsize=(4,2))
    ax = fig.add_subplot(111)
    im,ax = waterfall_plot(cube, ax)
    plt.colorbar(im, ax=ax)
    
    ax.set_xlabel('LSR Velocity (km/s)')
    ax.set_ylabel('Pixels')
    plt.show()
