# this code estimates the rms noise of a cube using different ways

from spectral_cube import SpectralCube
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from scipy.optimize import curve_fit
u.add_enabled_units(u.def_unit(['K (Tmb)'], represents=u.K))
u.add_enabled_units(u.def_unit(["K (Ta*)"], represents=u.K))

plt.rcParams['figure.dpi'] = 162

# class CubeNoise:
# def __init__(self, cube):
# self.cube = cube

def cube_trim(cube, vmin, vmax):
    cube = cube.spectral_slab(vmin*u.km/u.s, vmax*u.km/u.s)
    return cube

# vel_ranges of [-25,-15, 20,30] means -25 to -15 km/s AND 20 to 30 km/s
def rms_normal(cube,vel_ranges):
    # subcube = cube[:noise_channel,:,:] 
    # data = subcube.unmasked_data[:,:,:].value 
    slabs = []
    for idx in np.arange(0, len(vel_ranges), 2):
        
        slab = cube.spectral_slab(vel_ranges[idx] * u.km/u.s, vel_ranges[idx+1] * u.km/u.s)
        slabs.append(slab.unmasked_data[:,:,:].value)
    
    # Combine all selected slices into one array for calculation
    combined_data = np.concatenate(slabs, axis=0)    
    rms_map = np.nanstd(combined_data, axis=0)
    mean_rms = np.nanmean(rms_map)
    return rms_map, mean_rms

#for normal numpy cube, not spectral cube
# def rms_negative_np(npcube):
    # weighted_mean_cube_negative
    
def rms_negative(cube):
    if hasattr(cube, 'unmasked_data'):
        data = cube.unmasked_data[:,:,:].value
    else:
        data = cube
    masked_cube = np.where(data<0, data, np.nan)
    rms_map = np.sqrt(np.nanmean(masked_cube**2, axis=0))  
    mean_rms = np.nanmean(rms_map)
    max_rms = np.nanmax(rms_map)
    print(f'Maximum rms: {max_rms}')
    print(f'Mean rms: {mean_rms}')
    return rms_map, mean_rms
    
# this code has an issue where it doesn't work on data with nan values like CII    
def rms_histogram(cube):
    value = cube.unmasked_data[:,:,:].value.flatten()
    hist, bin_edges = np.histogram(value,
             bins=np.linspace(value.min(), value.max(), 200), density=True
             )
    bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
    sigma0 = np.std(value)
    mask = (bin_centers > -3*sigma0) & (bin_centers < 3*sigma0)

    p0 = [1.0, 0.0, sigma0]  # initial guess
    
    popt, pcov = curve_fit(gaussian, bin_centers[mask], hist[mask], p0=p0)
    return hist, bin_edges, bin_centers, popt
        
    
def gaussian(x, amp, mean, sigma):
    return amp * np.exp(-(x-mean)**2 / (2*sigma**2))
    
if __name__ == "__main__":
    cube_path = input("input full/relative path to cube fits: ")
    cube = SpectralCube.read(cube_path).with_spectral_unit(u.km/u.s)
    # cube_noise = CubeNoise(cube)
    print(f"Currently, the cube has a bandwidth of {cube.spectral_axis[0]:.2f} to {cube.spectral_axis[-1]:.2f} km/s")
    vlims = input("Input vmin and vmax in km/s separated by ',' to trim the cube (0 to skip trimming): ")
    if vlims != '0':
        vmin, vmax = list(map(float,vlims.split(',')))
        cube = cube_trim(cube, vmin, vmax)       
    option = int(input("select (1) for first N channel rms; (2) for negative channel rms; (3) for gaussian fitting of histogram: "))
    if option == 1:
        # channel = int(input("input number of channel to calculate rms: "))
        vel_ranges = input("input velocity ranges as emission free window: ")
        vel_ranges = list(map(float, vel_ranges.split(",")))
        rms_map, mean_rms = rms_normal(cube, vel_ranges)
        print(f'Mean rms noise: {mean_rms:.2f} K')
        plt.imshow(rms_map, origin='lower')
        plt.colorbar(label=cube.header['BUNIT'])
    elif option == 2:
        rms_map, mean_rms = rms_negative(cube)
        print(f'Mean rms noise: {mean_rms:.2f} K')
        # vmin = np.nanpercentile(rms_map,10)
        # vmax = np.nanpercentile(rms_map,100)
        percentile = float(input("Enter percentile to clip: "))
        vmax = np.nanpercentile(rms_map,percentile)
        # print(vmax)
        fig = plt.figure()  

        axis_format = input('format the axes in pixel or world coordinates?: ')
    
        if axis_format == 'world':
            wcs = cube.wcs.celestial
            ax = fig.add_subplot(111, projection=wcs)
        elif axis_format == 'pixel':
            ax = fig.add_subplot(111)
        else:
            print('unrecognised input, default to pixel coordinates')
            ax = fig.add_subplot(111)
        
        im = ax.imshow(rms_map, origin='lower', vmax=vmax)
        plt.colorbar(im, ax=ax, label=cube.header['BUNIT'])
    elif option == 3:
        hist, bin_edges, bin_centers, popt = rms_histogram(cube)
        amp, mean, sigma = popt
        print("Noise sigma =", sigma)
        plt.stairs(hist, bin_edges)
        # xfit = np.linspace(value.min(), value.max(), 200)
        xfit = bin_centers
        plt.plot(xfit, gaussian(xfit, *popt))
        plt.ylim(1e-5, 10)
        plt.axhline(amp/2,linestyle='dashed')
        plt.yscale('log')
        plt.ylabel('Count')
        plt.xlabel('T (K)')
        
<<<<<<< Updated upstream
    plt.show()
=======
    title = input('customise title? enter to skip: ')
    if title == '':
        ax.set_title('SNR Map')
    else:
        ax.set_title(title)
    plt.show()
    save_path = input('save the figure? type the path, or enter to skip: ')
    if save_path != '':
        fig.savefig(save_path, dpi=150, bbox_inches='tight', pad_inches=0)
        print(f'figure saved to {save_path}')
    else:
        print(f'figure not saved') 
>>>>>>> Stashed changes
