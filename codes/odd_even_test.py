# this codes is to test the average value of odd and even channels of each line-free pixel. we will use w43 c18o as an example.
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits




def odd_even_test(hdul):

    data = hdul['ORIGINAL'].data
    window = hdul['WINDOW'].data
    data_flat = data.reshape(data.shape[0],-1)
    ny, nx = data.shape[1], data.shape[2] 
    
    # a few methods to select spectra
    method = input('method: (1) white pixels (2) velocity range: ')

    if method == '1':
        mask = np.all(window==1, axis=0)
        pix_idx = np.flatnonzero(mask.ravel())
        line_free = data_flat[:, mask.ravel()]
    
    elif method == '2': 
        start = int(input('start channeL: '))
        end = int(input('end channeL: '))
        pix_idx = np.arange(nx * ny)
        line_free = data_flat[start:end]
    y, x = np.unravel_index(pix_idx, (ny, nx)) 
    odd_channels = line_free[::2,:]
    even_channels = line_free[1::2,:]
    
    odd_means = np.mean(odd_channels, axis=0)
    even_means = np.mean(even_channels, axis=0)

    return odd_means, even_means, x, y
    
    
if __name__ == '__main__':
    plt.rcParams.update({"font.size": 18,})
    fits_path = input("input path to the fits cube: ")
    hdul = fits.open(fits_path)
    odd_means,even_means,_,_ = odd_even_test(hdul)
    fig = plt.figure(tight_layout=True, figsize = (16,9))
    ax = fig.add_subplot(111)
    ax.scatter(odd_means,even_means)
    ax.set_xlabel('odd channels - mean')
    ax.set_ylabel('even channels - mean')
    ax.set_aspect('equal')
    ax.set_title(f"{hdul[0].header['OBJECT']} - {hdul[0].header['LINE']}")
    
    is_to_show = input('to show the plot? (y/n): ')
    if is_to_show == 'y':
        plt.show()
    
    output_path = input("input path to the output image (0 to skip): ")
    
    if output_path != '0':
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
