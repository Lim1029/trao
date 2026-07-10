# this code is to create a toy model to demonstrate the effect of "difference of gaussians" applied to a spectrum
# to access its potential in separating signal from background.

from codes.utils import spectrum_smooth
import astropy.units as u 
import matplotlib.pyplot as plt
from spectral_cube import SpectralCube
from matplotlib.widgets import Slider

# for now, have to manually define the cube path and the position  
cube_path = 'carta/TRAO/fits_files/w43/c18o.fits'
x,y = 52, 27

cube_path = 'carta/TRAO/fits_files/w43/c2h.fits'
x,y = 63, 26

cube = SpectralCube.read(cube_path).with_spectral_unit(u.km/u.s)

spectrum = cube[:,y,x].value
spectral_axis = cube.spectral_axis.value

# default parameters, will be changed in GUI slider later.
target_resolution_ori = 0.2
target_resolution_1 = 0.5
target_resolution_2 = 1.0
dof_multiplier = 1

# function to smooth a spectrum with two separate kernels and calculate their difference 
def calc_dof(target_resolution_1,target_resolution_2):
    smoothed_spectrum_1 = spectrum_smooth(spectrum, spectral_axis, target_resolution_1, resample=False)
    smoothed_spectrum_2 = spectrum_smooth(spectrum, spectral_axis, target_resolution_2, resample=False)
    
    difference_of_gaussians = smoothed_spectrum_1 - smoothed_spectrum_2

    return smoothed_spectrum_1, smoothed_spectrum_2, difference_of_gaussians

# ================= plot =================
fig, ax = plt.subplots(figsize=(8, 5))
plt.subplots_adjust(left=0.25, bottom=0.45)
smoothed_spectrum_ori = spectrum_smooth(spectrum, spectral_axis, target_resolution_ori, resample=False)
smoothed_spectrum_1, smoothed_spectrum_2, difference_of_gaussians = calc_dof(target_resolution_1, target_resolution_2)
(line_ori,) = ax.step(spectral_axis, smoothed_spectrum_ori, label='Ori')
(line_g1,) = ax.step(spectral_axis, smoothed_spectrum_1, label='G1')
(line_g2,) = ax.step(spectral_axis, smoothed_spectrum_2, label='G2')
(line_dof,) = ax.step(spectral_axis, difference_of_gaussians*dof_multiplier, label='DoG')
ax.set_xlabel("Velocity (km/s)")
ax.set_ylabel("TA* (K)")
ax.legend()

# ================= sliders =================
axcolor = "lightgoldenrodyellow"

def add_slider(ypos, label, vmin, vmax, vinit):
    axs = plt.axes([0.25, ypos, 0.65, 0.03], facecolor=axcolor)
    return Slider(axs, label, vmin, vmax, valinit=vinit)

s_smori= add_slider(0.34, "Ori Spec smooth (km/s)", 0.1, 5, 0.2)
s_gk1 = add_slider(0.30, "G1 σ (km/s)", 0.1, 30, 0.5)
s_gk2 = add_slider(0.26, "G2 σ (km/s)", 0.1, 30, 1.0)
s_dofmul = add_slider(0.22, "dof multiplier", 1, 20, 1.0)

# ================= update function =================
# update the smoothed spectra and dof-ed spectrum
def update(val):    
    smoothed_spectrum_1, smoothed_spectrum_2, difference_of_gaussians = calc_dof(s_gk1.val, s_gk2.val)
    line_g1.set_ydata(smoothed_spectrum_1)
    line_g2.set_ydata(smoothed_spectrum_2)
    line_dof.set_ydata(difference_of_gaussians*s_dofmul.val)
    # ax.relim()            # recompute data limits, not necessary for now
    # ax.autoscale_view()   # rescale the plot to the new data limits, not necessary for now
    fig.canvas.draw_idle()

# update the original reference spectrum
def update_ori(val):
    smoothed_spectrum_ori = spectrum_smooth(spectrum, spectral_axis, s_smori.val, resample=False)
    line_ori.set_ydata(smoothed_spectrum_ori)

for s in [s_gk1, s_gk2, s_dofmul]:
    s.on_changed(update)

s_smori.on_changed(update_ori)

plt.show()
