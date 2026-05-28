# Baseline Method

Here is the discription about the method we are using for the baseline subtraction.

The overall workflow for the baseline subtration can be understood by this flowchart

<img width="600" height="900" alt="base_flow" src="https://github.com/user-attachments/assets/192e80ce-16ee-4e6c-b372-3e3c63a2413e" />

## Preprocessing the Cube

Optionally the the spectral cube can be trimmed to a velocity range to remove the bad channels in either end of the spectra.

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
