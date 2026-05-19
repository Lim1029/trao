from astropy.convolution import convolve_fft
import numpy as np
from scipy.ndimage import uniform_filter
# adopted from BTS code 
# https://github.com/SeamusClarke/BTS
# different from BTS is that the Filter_size is now a user defined set (v,y,x)
def TopHat_3DFilter(Image, Filter_size):
    kernel = np.ones(Filter_size)
    kernel = kernel / np.sum(kernel)
    # Final_image = convolve_fft(Image,kernel,boundary="wrap")
    # alternative
    Final_image = uniform_filter(Image, size=Filter_size, mode='wrap')
    return Final_image
