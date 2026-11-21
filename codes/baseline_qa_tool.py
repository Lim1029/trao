#  the tool should be able to visualise window, original, baselined, and the model, and the ability to hide/show any. also need to make moment-0 map of original and baselined cube, of all channels and window, non-window channel. I do not know if this can be scripted in CARTA then no need to built the plot. should also have a comment textbox to show whatever comments the cube might have (at primary hdu)... also need to have an option to show average spectrum (for now, across all cube is okay, localise averaging would be complicated) i should code myself first then get ai to help.

fits_path = input("Input full/relative path to cube fits: ")
from astropy.io import fits
from astropy.units import Unit
import astropy.units as u 
import numpy as np
from utils.spectrum_smooth import spectrum_smooth_kernel
from matplotlib.widgets import Button, TextBox
import matplotlib.pyplot as plt
import textwrap
hdul = fits.open(fits_path)
hdu_names = ['ORIGINAL','BASELINED','WINDOW','BASELINE']

plt.rcParams.update({"font.size": 8,})

# fig = plt.figure(tight_layout=True, figsize = (16,9))
fig = plt.figure(figsize = (16,6))
gs = fig.add_gridspec(5,1, height_ratios=[0,1.5,0.25,0.1,0.2])
gs_mom0 = gs[0, 0].subgridspec(1, 6)

# we assume the velocity axis of each hdu is the same
# so we derive the axes from their respective header
header = hdul[0].header
pix_axis = np.arange(header['NAXIS3'])
vel = header['CRVAL3'] + (pix_axis - (header['CRPIX3'] - 1.0)) * header['CDELT3']
vel = vel*Unit(header['CUNIT3']).to(u.km/u.s)

data = {}
for hdu_name in hdu_names:
    if hdu_name == 'BASELINED':
        data['BASELINED'] = hdul[0].data
    else:
        data[hdu_name] = hdul[hdu_name].data

# make moment0
mom0_original_all = np.mean(data['ORIGINAL'],axis=0)
mom0_original_window = np.nanmean(np.where(data['WINDOW']==0,data['ORIGINAL'],np.nan),axis=0)
mom0_original_nonwindow = np.nanmean(np.where(data['WINDOW']==1,data['ORIGINAL'],np.nan),axis=0)
mom0_baselined_all = np.mean(data['BASELINED'],axis=0)
mom0_baselined_window = np.nanmean(np.where(data['WINDOW']==0,data['BASELINED'],np.nan),axis=0)
mom0_baselined_nonwindow = np.nanmean(np.where(data['WINDOW']==1,data['BASELINED'],np.nan),axis=0)

# fig for moment0
fig_mom = plt.figure(figsize = (16,3))

ax_mom0_original_all = fig_mom.add_subplot(161)
ax_mom0_original_window = fig_mom.add_subplot(162)
ax_mom0_original_nonwindow = fig_mom.add_subplot(163)
ax_mom0_baselined_all = fig_mom.add_subplot(164)
ax_mom0_baselined_window = fig_mom.add_subplot(165)
ax_mom0_baselined_nonwindow = fig_mom.add_subplot(166)

im = ax_mom0_original_all.imshow(mom0_original_all, origin='lower')
im.set_rasterized(True)
plt.colorbar(im, ax=ax_mom0_original_all, location='top')
im = ax_mom0_original_window.imshow(mom0_original_window, origin='lower')
im.set_rasterized(True)
plt.colorbar(im, ax=ax_mom0_original_window, location='top')
im = ax_mom0_original_nonwindow.imshow(mom0_original_nonwindow, origin='lower')
im.set_rasterized(True)
plt.colorbar(im, ax=ax_mom0_original_nonwindow, location='top')
im = ax_mom0_baselined_all.imshow(mom0_baselined_all, origin='lower')
im.set_rasterized(True)
plt.colorbar(im, ax=ax_mom0_baselined_all, location='top')
im = ax_mom0_baselined_window.imshow(mom0_baselined_window, origin='lower')
im.set_rasterized(True)
plt.colorbar(im, ax=ax_mom0_baselined_window, location='top')
im = ax_mom0_baselined_nonwindow.imshow(mom0_baselined_nonwindow, origin='lower')
im.set_rasterized(True)
plt.colorbar(im, ax=ax_mom0_baselined_nonwindow, location='top')

# write titles
ax_mom0_original_all.set_title('Original; all',pad=30,fontsize=8)
ax_mom0_original_window.set_title('Original; window only',pad=30,fontsize=8)
ax_mom0_original_nonwindow.set_title('Original; non-window only',pad=30,fontsize=8)
ax_mom0_baselined_all.set_title('Baselined; all',pad=30,fontsize=8)
ax_mom0_baselined_window.set_title('Baselined; window only',pad=30,fontsize=8)
ax_mom0_baselined_nonwindow.set_title('Baselined; non-window only',pad=30,fontsize=8)

# plot markers placeholders
px_marker_1, = ax_mom0_original_all.plot([], [], marker='x', color='red', linestyle='none', ms=8, mew=1.5)
px_marker_2, = ax_mom0_original_window.plot([], [], marker='x', color='red', linestyle='none', ms=8, mew=1.5)
px_marker_3, = ax_mom0_original_nonwindow.plot([], [], marker='x', color='red', linestyle='none', ms=8, mew=1.5)
px_marker_4, = ax_mom0_baselined_all.plot([], [], marker='x', color='red', linestyle='none', ms=8, mew=1.5)
px_marker_5, = ax_mom0_baselined_window.plot([], [], marker='x', color='red', linestyle='none', ms=8, mew=1.5)
px_marker_6, = ax_mom0_baselined_nonwindow.plot([], [], marker='x', color='red', linestyle='none', ms=8, mew=1.5)


def update_mom_marker(x,y):
    px_marker_1.set_data([x], [y])
    px_marker_2.set_data([x], [y])
    px_marker_3.set_data([x], [y])
    px_marker_4.set_data([x], [y])
    px_marker_5.set_data([x], [y])
    px_marker_6.set_data([x], [y])
    fig_mom.canvas.draw_idle() 


def on_click(event):
    if event.inaxes not in [ax_mom0_original_all, ax_mom0_original_window, ax_mom0_original_nonwindow, ax_mom0_baselined_all, ax_mom0_baselined_window, ax_mom0_baselined_nonwindow]:
        return
    if event.xdata is None or event.ydata is None:
        return
    # ignore clicks while zoom/pan tools are active
    tb = fig_mom.canvas.toolbar
    if tb is not None and tb.mode != '':
        return

    x_ = int(np.clip(np.round(event.xdata), 0, nx - 1))
    y_ = int(np.clip(np.round(event.ydata), 0, ny - 1))

    state['x'], state['y'] = x_, y_
    update_plot(x_, y_)

fig_mom.canvas.mpl_connect('button_press_event', on_click)

# perform odd-even test
from odd_even_test import odd_even_test
odd_means,even_means, x_oet, y_oet = odd_even_test(hdul)

fig_oet = plt.figure(figsize=(5,5))
ax_oet = fig_oet.add_subplot(1,1,1)
ax_oet.set_aspect("equal")
sc = ax_oet.scatter(odd_means, even_means, s=1, picker=True)
ax_oet.set_xlabel('odd mean')
ax_oet.set_ylabel('even mean')

def on_pick(event):
    if event.artist is not sc:
        return

    # several points can fall within the pick radius; take the closest one
    ind = event.ind
    mx, my = event.mouseevent.xdata, event.mouseevent.ydata
    dists = (odd_means[ind] - mx)**2 + (even_means[ind] - my)**2
    i = ind[np.argmin(dists)]
    px, py = x_oet[i], y_oet[i]
    state['x'], state['y'] = px, py

    update_plot(px, py)

cid = fig_oet.canvas.mpl_connect('pick_event', on_pick)


# plot the average spectrum
avg_original = np.mean(data['ORIGINAL'],axis=(1,2))
avg_baselined = np.mean(data['BASELINED'],axis=(1,2))

# plot the average (median) white-pixels spectrum. white-pixel here means pixels with no emission window
mask = np.all(data['WINDOW'] == 1, axis=0) 
spec_flat = data['ORIGINAL'].reshape(data['ORIGINAL'].shape[0],-1)
selected = spec_flat[:, mask.ravel()]
avg_white_original = np.median(selected, axis=1)
std_white = np.std(avg_white_original)

x,y = 33,33
ax = fig.add_subplot(gs[1,0])
ax.autoscale(enable=True, axis='y')
v0 = vel[0]
dv = vel[1] - vel[0]

def vel2chan(v):
    return (v - v0) / dv

def chan2vel(c):
    return v0 + c * dv
ax_top = ax.secondary_xaxis('top', functions=(vel2chan, chan2vel))
_l_avg_original, = ax.step(vel, avg_original, linewidth=0.5, color='pink', label='AVG ORIGINAL')
_l_avg_baselined, = ax.step(vel, avg_baselined, linewidth=0.5, color='cyan', label='AVG BASELINED')
_l_avg_white, = ax.step(vel, avg_white_original, linewidth=0.5, color='black', label='AVG WHITE')

# rank the pixels based on original all mom0
flat_max = mom0_original_all.ravel()
order    = np.argsort(flat_max)   
order = order[::-1] # strongest first
iy_all, ix_all = np.unravel_index(order, mom0_original_all.shape)
rank = np.where((iy_all == y) & (ix_all == x))[0][0] + 1
rank_map = np.empty(mom0_original_all.shape, dtype=int)
rank_map[iy_all, ix_all] = np.arange(1, mom0_original_all.size + 1)
total = mom0_original_all.size

def setup_plot():
    global _l_window, _rank_text, _l_orig, _l_base, _l_baselined, _spans, _stats_text

    # rank / coordinate annotation
    _rank_text = ax.text(0, 1, '', transform=ax.transAxes,
                         ha='left', va='bottom')

    _stats_text = ax.text(0.05, 0.95, '', transform=ax.transAxes,
                         ha='left', va='top')

    # three step-spectra lines (kept as Line2D for repeated set_data)
    _l_orig,       = ax.step([], [], lw=0.5, label='ORIGINAL')
    _l_base,       = ax.step([], [], lw=1.0, label='BASELINE')
    _l_baselined,  = ax.step([], [], lw=0.5, label='BASELINED')
    _l_window,     = ax.plot([], [], lw=0.5, label='WINDOW',)



# add a row for showing comments
comments = hdul[0].header['COMMENT']
ax_comment = fig.add_subplot(gs[4,0])
comment_string = ''.join(comments)
wrapped = textwrap.fill(comment_string, width=100)
ax_comment.text(0, 1, wrapped, ha='left', va='top', transform=ax_comment.transAxes)
ax_comment.axis('off')

gs_control = gs[2,0].subgridspec(1,3,width_ratios=[1,0.5,1])

# add pixel navigator
gs_navi = gs_control[0,0].subgridspec(1,4)
ax_left = fig.add_subplot(gs_navi[0,0])
ax_right = fig.add_subplot(gs_navi[0,1])
ax_up = fig.add_subplot(gs_navi[0,2])
ax_down = fig.add_subplot(gs_navi[0,3])

# add toggle navigator
gs_toggle = gs[3,0].subgridspec(1,7)
ax_original_toggle = fig.add_subplot(gs_toggle[0,0])
ax_baselined_toggle = fig.add_subplot(gs_toggle[0,1])
ax_avg_original_toggle = fig.add_subplot(gs_toggle[0,2])
ax_avg_baselined_toggle = fig.add_subplot(gs_toggle[0,3])
ax_avg_white_toggle = fig.add_subplot(gs_toggle[0,4])
ax_baseline_toggle = fig.add_subplot(gs_toggle[0,5])
ax_relim = fig.add_subplot(gs_toggle[0,6])

# add smooth textbox
gs_smooth = gs_control[0,2].subgridspec(3,2) 
ax_smooth_original = fig.add_subplot(gs_smooth[0,0])
ax_smooth_baselined = fig.add_subplot(gs_smooth[1,0])
ax_smooth_white = fig.add_subplot(gs_smooth[2,0])
ax_smooth = fig.add_subplot(gs_smooth[:,1])

class FastTextBox(TextBox):
    def _rendercursor(self):
        drawon, self.drawon = self.drawon, False   # suppress the blocking draw()
        super()._rendercursor()
        self.drawon = drawon
        self.ax.figure.canvas.draw_idle()          # coalesced redraw instead

kw = dict(initial='10', color='0.95', hovercolor='0.95')
tb_smooth_original  = FastTextBox(ax_smooth_original,  'kernel size [original] (pix)',  **kw)
tb_smooth_baselined = FastTextBox(ax_smooth_baselined, 'kernel size [baselined] (pix)', **kw)
tb_smooth_white     = FastTextBox(ax_smooth_white,     'kernel size [white] (pix)',     **kw)

# tb_smooth_original = TextBox(ax_smooth_original, 'kernel size [original] (pix)', initial='0')
# tb_smooth_baselined = TextBox(ax_smooth_baselined, 'kernel size [baselined] (pix)', initial='0')
# tb_smooth_white = TextBox(ax_smooth_white, 'kernel size [white] (pix)', initial='0')
b_smooth = Button(ax_smooth, "smooth")

b_left = Button(ax_left, '◀')
b_right = Button(ax_right, '▶')
b_up = Button(ax_up, '▲')
b_down = Button(ax_down, '▼')

b_original_toggle = Button(ax_original_toggle, 'Show/Hide Original')
b_baselined_toggle = Button(ax_baselined_toggle, 'Show/Hide Baselined')
b_avg_original_toggle = Button(ax_avg_original_toggle, 'Show/Hide Original Average')
b_avg_baselined_toggle = Button(ax_avg_baselined_toggle, 'Show/Hide Baselined Average')
b_avg_white_toggle = Button(ax_avg_white_toggle, 'Show/Hide White Average')
b_baseline_toggle = Button(ax_baseline_toggle, 'Show/Hide Baseline Model')

b_relim = Button(ax_relim, 'axis relim')

def toggle_original(val):
    _l_orig.set_visible(not _l_orig.get_visible())
    fig.canvas.draw_idle()

def toggle_baselined(val):
    _l_baselined.set_visible(not _l_baselined.get_visible())
    fig.canvas.draw_idle()

def toggle_avg_original(val):
    _l_avg_original.set_visible(not _l_avg_original.get_visible())
    fig.canvas.draw_idle()

def toggle_avg_baselined(val):
    _l_avg_baselined.set_visible(not _l_avg_baselined.get_visible())
    fig.canvas.draw_idle()

def toggle_avg_white(val):
    _l_avg_white.set_visible(not _l_avg_white.get_visible())
    fig.canvas.draw_idle()

def toggle_baseline(val):
    _l_base.set_visible(not _l_base.get_visible())
    fig.canvas.draw_idle()

def relim(val):
    ax.relim()            
    ax.autoscale_view() 
    fig.canvas.draw_idle()

b_original_toggle.on_clicked(toggle_original)
b_baselined_toggle.on_clicked(toggle_baselined)
b_avg_original_toggle.on_clicked(toggle_avg_original)
b_avg_baselined_toggle.on_clicked(toggle_avg_baselined)
b_baseline_toggle.on_clicked(toggle_baseline)
b_avg_white_toggle.on_clicked(toggle_avg_white)
b_relim.on_clicked(relim)


state = {'x': x, 'y': y}
ny, nx = data['ORIGINAL'].shape[1], data['ORIGINAL'].shape[2]

def update_spectrum(x,y):
    # 2) coordinate + rank text
    rank = rank_map[y, x]
    _rank_text.set_text(f'(x,y)=({x},{y}); mom0 rank = ({rank}/{total})')

    # calculate stats
    original = data['ORIGINAL'][:,   y, x]
    mask = data['WINDOW'][:, y, x]
    linefree = np.where(mask==1, original, np.nan)
    linebusy = np.where(mask==0, original, np.nan)
    std_linefree = np.nanstd(linefree)
    std_linebusy = np.nanstd(linebusy)
    std_all = np.nanstd(original)
    _stats_text.set_text(f'Std Original (all, emission, non-emission, white) = ({std_all:.2e}, {std_linebusy:.2e}, {std_linefree:.2e}, {std_white:.2e})')
    
    # smooth the spectra first if necessary:
    if tb_smooth_original.text != '0':
        kernel_size = int(tb_smooth_original.text)
        smoothed_original = spectrum_smooth_kernel(data['ORIGINAL'][:,   y, x], vel, kernel_size)
        _l_orig.set_data(vel, smoothed_original)
    else:
        _l_orig.set_data(vel,   data['ORIGINAL'][:,   y, x])

    if tb_smooth_baselined.text != '0':
        kernel_size = int(tb_smooth_baselined.text)
        smoothed_baselined = spectrum_smooth_kernel(data['BASELINED'][:,   y, x], vel, kernel_size)
        _l_baselined.set_data(vel, smoothed_baselined)
    else:
        _l_baselined.set_data(vel,  data['BASELINED'][:,  y, x])

    if tb_smooth_white.text != '0':
        kernel_size = int(tb_smooth_white.text)
        smoothed_white = spectrum_smooth_kernel(avg_white_original, vel, kernel_size)
        _l_avg_white.set_data(vel, smoothed_white)
    else:
        _l_avg_white.set_data(vel, avg_white_original)

    _l_base.set_data(     vel,   data['BASELINE'][:,   y, x])
    _l_window.set_data(   vel,   1-data['WINDOW'][:,   y, x])

    fig.canvas.draw_idle()


# ===== per-pixel update =====
def update_plot(x, y):
    update_mom_marker(x, y)
    update_spectrum(x,y)


def move(dx, dy):
    x_ = int(np.clip(state['x'] + dx, 0, nx - 1))
    y_ = int(np.clip(state['y'] + dy, 0, ny - 1))
    if (x_, y_) != (state['x'], state['y']):
        #print('position state changed')
        state['x'], state['y'] = x_, y_
        update_plot(x_, y_)

# after creating all axes:
setup_plot()
ax.legend()
update_plot(state['x'], state['y'])

b_up   .on_clicked(lambda e: move(0,  1))
b_down .on_clicked(lambda e: move(0, -1))
b_left .on_clicked(lambda e: move(-1, 0))
b_right.on_clicked(lambda e: move( 1, 0))

b_smooth.on_clicked(lambda e: update_spectrum(state['x'], state['y']))


fig.suptitle(f"{hdul[0].header['OBJECT']} - {hdul[0].header['LINE']}", fontsize=24)

plt.show()
