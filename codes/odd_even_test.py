# this codes is to test the average value of odd and even channels of each line-free pixel. we will use w43 c18o as an example.
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

fits_path = input("input path to the fits cube: ")

plt.rcParams.update({"font.size": 18,})

hdul = fits.open(fits_path)

fig = plt.figure(tight_layout=True, figsize = (16,9))

data = hdul['ORIGINAL'].data
window = hdul['WINDOW'].data

# mask the data to only line-free pixels

mask = np.all(window==1, axis=0)
data_flat = data.reshape(data.shape[0],-1)
line_free = data_flat[:, mask.ravel()]
odd_channels = line_free[::2,:]
even_channels = line_free[1::2,:]

# odd_channels = data_flat[::2,:]
# even_channels = data_flat[1::2,:]
odd_means = np.mean(odd_channels, axis=0)
even_means = np.mean(even_channels, axis=0)

# plt.plot(np.nanmean(data_flat, axis=1))

plt.scatter(odd_means,even_means)
plt.xlabel('odd channels - mean')
plt.ylabel('even channels - mean')
plt.title(f"{hdul[0].header['OBJECT']} - {hdul[0].header['LINE']}")
plt.show()
