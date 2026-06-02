# Baseline Method


The cube_baseline_subtraction is being run as a module here from the directory which contains codes/cube_baseline_subtraction.py

Here is the description about the method we are using for the baseline subtraction.

The overall workflow for the baseline subtraction can be understood by this flowchart

<img width="600" height="900" alt="base_flow" src="https://github.com/user-attachments/assets/192e80ce-16ee-4e6c-b372-3e3c63a2413e" />

##  # Preprocessing the Cube

### Trimming 

Optionally the the spectral cube can be trimmed to a velocity range to remove the bad channels in either end of the spectra.

Parameters : 

    vmin ( minimum velocity limit)
    vmax ( upper velocity limit)

### Smoothing

And also the spectra can be smoothed to desired resolution to get better snr for signal detection.

Parameters :

      target_reso ( Resolution to which to smooth the spectra)


##  # Masking Methods

There are three option in present pipeline for masking 

    1. No masking
    2. Fixed Window
    3. Automatic windowing

### No window

In this case we are defining that all the channels for a pixel is baseline, there is no emission part.

### Fixed window

We define a velocity range to consider as emission. 
for example if the spectral range is from 0 to 160 km/s and we are defining that consider 50 to 60 km/s as emission, then the channels outside 50-60 km/s will be assigned as 1 (baseline).

parameters:

    vlims (pair of velocity ranges to mask as emission)

### Automatic window

One other masking method we have is automatic windowing adapted from [this paper](https://ui.adsabs.harvard.edu/abs/2015MNRAS.453...73J/abstract).

The python code is in the file 

    codes/cube_signal_window.py

Additionally as a option for weak signal we can also smooth the spectra with a 3D smoothing with smooth_kernel to get better SNR. for this we are using the 

    codes/utils/cube_smooth_tophat.py 

The working of this method can be understood by this [illustration](https://lim1029.github.io/jcmt_window_demo.pdf)

parameters:

    nbin ( number of bins to divide the channels )
    clips ( sigma levels to estimate the threshold to identify the emission bins (currently we are using three iterations)
    smooth_kernel ( used for 3D smoothing the cube ) 
    bin_expand ( number of neighbouring bins to expand to consider as emission)
    

## # Baseline Fitting methods

So for the baseline fitting we also have three methods for now 

    1. Fixed polynomial
    2. Iterative Polynomial
    3. Spline Fit

### Fixed Polynomial

In this method for every pixel a single fixed order(user defined) polynomial is used to fit the baseline in the part defined as baseline by the mask created in the above step.

Then the baseline model is subtracted from the raw spectra for each pixel.

parameters: 

    poly_order ( polynomial order to fit baseline )

### Iterative Polynomial 

For this method for each pixel we take the following steps
    
    we start with poly order = 0
    fit the baseline and compute AIC (Akaike information criterion ) which indicate the goodness of fit
    then increase the order by 1 and compute AIC again 
    the iteration stop when the AIC starts to get worse.
    the order with lowest AIC is chosen to fit the baseline and subtracted from the spectra 

An upper limit to the polynomial order is required to define first.

parameters :

    max_order ( maximum order for iteration )

For AIC this note can be referred : [AIC_statistics.pdf](https://github.com/user-attachments/files/28490685/AIC_statistics.pdf)


### Spline Fit

By default order 3 ( cubic spline ) is used.

parameters

    knot_start ( channel index from where the interior knot begins)
    knot_spacing ( channel spacing between 2 knots )
    edge_channels ( number of edge channels to take a mean or median to get a straight line to replace the emission part. )

For an illustration of how these parameters affect the spline fitting, refer to [spline_illustration](https://github.com/vinay-ydv19/w43-doc/blob/main/Spline_Baseline.pdf)

short illustration to understand how individual splines are constructed between knots, and summed to produce the baseline [Making sense of B spline.pdf](https://github.com/user-attachments/files/28490718/Making.sense.of.B.spline.pdf)


## # Output FITS File

Primary HDU:
    Baseline-subtracted cube

BASELINE:
    Fitted baseline model

UN-BASELINED:
    Original cube before subtraction

MASK:
    Spectral mask used during fitting



## Optional Visualization

After baseline subtraction the average spectra can be plotted to visualise the result. 

plotting code 

    codes/utils/average_plotting.py

we have two option for plotting either all the pixels can be averaged or pixels with above certain SNR threshold can be averaged.

example of average plot
<img width="800" height="450" alt="w40_c18o_spline" src="https://github.com/user-attachments/assets/b9a7d2dc-bd3b-411f-bcf5-b3db357a70a3" />


## # How to Run 

To perform baseline subtraction this terminal command can be used 

    python -m codes.cube_baseline_subtraction


## Example

    Input full/relative path to cube fits: 
            /home/vinay/narit/w43/fits/trao/w43_c18o_mos.fits

    Input vmin and vmax in km/s separated by ',' to trim the cube (0 to skip trimming):
            20,150
    
    Input target resolution (km/s) to perform smoothing prior to baseline fitting (0 to skip):
            0

    How to define spectral window? (1) no window (2) user-defined (3) automatic: 
            3

    Input nbin, to divide 3121 channels: 
            40

    Input clip (in unit of sigma) separated by ',':  
            1.5,2,2.5
    
    tophat smooth cube prior to windowing? (enter kernel size v,y,x or 0 to skip): 
            5,5,5

    Number of neighboring bin each side to mask as emission to protect wings: 
            2

    Input Baseline Fitting method: (1) fixed poly (2) iterative poly (3) spline: 
            3

    Input interior knot starting channel: 
            200

    Input knot spacing in channels: 
            200

    Input edge averaging channels: 
            50

    Visualise baseline fitting? (y/n): 
            y

    Enter SNR threshold for cube masking (0 to skip): 
            0

    Save this plot? (y/n): 
            n

    input the output path: 
            /home/vinay/narit/w43/fits/trao/base/spline/c18o.fits
    
    
    


