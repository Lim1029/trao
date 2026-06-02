import numpy as np

# this function is obtained from https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.BSpline.html
def B(x, k, i, t):
	if k == 0:
		return 1.0 if t[i] <= x < t[i+1] else 0.0
	if t[i+k] == t[i]:
		c1 = 0.0
	else:
		c1 = (x - t[i])/(t[i+k] - t[i]) * B(x, k-1, i, t)
	if t[i+k+1] == t[i+1]:
		c2 = 0.0
	else:
		c2 = (t[i+k+1] - x)/(t[i+k+1] - t[i+1]) * B(x, k-1, i+1, t)
	return c1 + c2

def bspline(x, t, c, k):
	n = len(t) - k - 1
	assert (n >= k+1) and (len(c) >= n)
	return sum(c[i] * B(x, k, i, t) for i in range(n))

import matplotlib.pyplot as plt
plt.style.use('./codes/astro.mplstyle')
plt.rcParams['figure.figsize'] = (10,5)
x_list = np.arange(0,18,1) 

# here we modified the order of the spline (k) to generate
k = 2

# we hardocde the knot positions to be equidistance along the x-axis
t = [0,3,6,9,12,15,18]

# we hardcode the knot positions, which is 4 knots at both ends, and equidistance in between.
# t = [0,0,0,0,3,6,9,12,15,18,18,18,18]

n = len(t)-k-1 # the number of spline

# here we can choose whether to set the spline coefficient as the same, or random
c = [1 for i in range(n)] # same spline coefficient
# random spline coefficient c
# rng = np.random.default_rng()
# c = rng.integers(low=-5, high=6, size=n)

spline_list = []
for i in range(n):
	spline = []
	for x in x_list:
		spline.append(B(x,k,i,t))
	spline_list.append(spline)

spline_list = np.array(spline_list)
count = 0
for i,spline in enumerate(spline_list):
	spline = np.array(spline)
	plt.plot(x_list,spline*c[i], label='spline at a point')
	count += 0
c = np.repeat(c, len(x_list)).reshape(spline_list.shape)
total_spline = np.sum(spline_list*c, axis=0)
plt.plot(x_list, total_spline, color='black', label='total spline')
plt.scatter(t,np.zeros_like(t),label='knots')
plt.legend()
plt.show()
