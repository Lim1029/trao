# Baseline Method

Here is the discription about the method we are using for the baseline subtraction.

The overall workflow for the baseline subtration can be understood by this flowchart

<img width="600" height="900" alt="base_flow" src="https://github.com/user-attachments/assets/192e80ce-16ee-4e6c-b372-3e3c63a2413e" />

## Preprocessing the Cube

### Trimming 

Optionally the the spectral cube can be trimmed to a velocity range to remove the bad channels in either end of the spectra.

Parameters : vmin ( minimum velcoity limit)
             vmax ( upper velocity limit)

And also the spectra can be smoothed to desired resolution to get better snr for signal ditection.


## Masking Methods

There are three option in present pipeline for masking 

    1. No masking
    2. Fixed Window
    3. Automatic windowing

### No window

In this case we are defining that all the channels for a pixel is baseline, there is no emission part.

### Fixed window

We define a velocity range to consider as emission. 
for example if the spectral range is from 0 to 160 km/s and we are defining that consider 50 to 60 km/s as emission, then the channels outside 50-60 km/s will be assigned as 1 (baseline).

### Automatic window

One other masking method we have is automatic windowing adapted from [this paper](https://ui.adsabs.harvard.edu/abs/2015MNRAS.453...73J/abstract).

The working of this method can be understood by this [illustration](https://lim1029.github.io/jcmt_window_demo.pdf)

## Baseline Fitting methods

So for the baseline fitting we also have three methods for now 

    1. Fixed polynomial
    2. Iterative Polynomial
    3. Spline Fit

### Fixed Polynomial

In this method for every pixel a single fixed order(user defined) polynomial is used to fit the baseline in the part defined as baseline by the mask created in the above step.

Then the baseline model is subtracted from the raw spectra for each pixel.

### Iterative Polynomial 

For this method for each pixel we take the following steps
    
    we start with poly order = 0
    fit the baseline and compute AIC (Akaike information criterion ) which indicate the goodness of fit
    then increase the order by 1 and compute AIC again 
    the iteration stop when the AIC starts to get worse.
    the order with lowest AIC is choosen to fit the baseline and subtracted from the spectra 


### Spline Fit


