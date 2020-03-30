linear_wave_test = True
sod_shock_test = False

import matplotlib
#matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib import colors as colors
import matplotlib.gridspec as gridspec
from matplotlib.ticker import FormatStrFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import rc
from matplotlib import patches
import matplotlib.lines as mlines
#rc("text", usetex=True)
rc("font", family="serif")
rc("font", size=14)
rc("axes", grid=True)
rc("grid", linestyle="--")
rc("xtick", direction="in")
rc("ytick", direction="in")
rc("savefig", format="pdf", bbox="tight")
import numpy as np
import os
import multiprocessing

param = np.loadtxt("data/param.dat")
nx = param[0]
ny = param[1]
dx = param[2]
dy = param[3]
xmin = param[4]
xmax = param[5]
ymin = param[6]
ymax = param[7]
nfile = int(param[8])
img_rat = (xmax - xmin) / (ymax - ymin);

m_air =	1.204e-3 / 0.02504e21
kB = 1.38e-16

def doit(fnum):
	print("fnum = %d" % fnum)
	try:
		rho = np.transpose(np.loadtxt("data/rho_%05d.dat" % fnum))
		vx = np.transpose(np.loadtxt("data/vx_%05d.dat" % fnum))
		vy = np.transpose(np.loadtxt("data/vy_%05d.dat" % fnum))
		press = np.transpose(np.loadtxt("data/press_%05d.dat" % fnum))
		v = np.sqrt(vx**2 + vy**2)
	except:
		return

	x = np.linspace(xmin, xmax, nx+1)
	xcc = (x[:-1] + x[1:]) / 2
	y = np.linspace(ymin, ymax, ny+1)
	ycc = (y[:-1] + y[1:]) / 2
	[x, y] = np.meshgrid(x, y, indexing="xy")

	def pl(ax, q, cmap="viridis", vbound=None):
		if vbound is None:
			vmin = 0
			vmax = np.amax(q)
		else:
			[vmin, vmax] = vbound

		im = ax.pcolormesh(x, y, q, cmap=cmap, norm=colors.Normalize(vmin=vmin, vmax=vmax))
		ax.set_aspect(1)
		divider = make_axes_locatable(ax)
		cax = divider.append_axes("right", size="5%", pad=0.05)
		cbar = plt.colorbar(im, cax=cax)
		cbar.ax.tick_params()

	if linear_wave_test:
		fig = plt.figure(figsize=(15, 7.5), dpi=72)
		gs = gridspec.GridSpec(1, 2)

		ax = fig.add_subplot(gs[0,0])
		pl(ax, rho - 1, cmap="cividis", vbound=[-1e-6, 1e-6])

		ax = fig.add_subplot(gs[0,1])
		pl(ax, press - 3.0/5.0, cmap="cividis", vbound=[-1e-6, 1e-6])
	elif sod_shock_test:
		fig = plt.figure(figsize=(15, 15), dpi=72)
		gs = gridspec.GridSpec(2, 2)

		rho_cmap = "GnBu"
		ax = fig.add_subplot(gs[0,0])
		pl(ax, rho, cmap=rho_cmap, vbound=[0,2.25])

		t_cmap = "YlOrBr_r"
		ax = fig.add_subplot(gs[0,1])
		pl(ax, press / rho, cmap=t_cmap, vbound=[0,1.5])

		ax = fig.add_subplot(gs[1,0])
		ax.plot(xcc, rho[int(ny/2),:], "C0.-")

		ax = fig.add_subplot(gs[1,1])
		ax.plot(xcc, vx[int(ny/2),:], "C1.-")
	else:
		fig = plt.figure(figsize=(15, 7.5), dpi=72)
		gs = gridspec.GridSpec(1, 2)

		rho_cmap = "GnBu"
		ax = fig.add_subplot(gs[0,0])
		#pl(ax, np.log10(rho), cmap=rho_cmap, vbound=[-8,1])
		pl(ax, rho, cmap=rho_cmap, vbound=[0,2.25])
		#ax.streamplot(xcc, ycc, vx, vy, color=v, cmap="gray", density=2)
		#res_circle = plt.Circle((0, 0), 1.2e10, color="w", fill=False, linewidth=0.5)
		#ax.add_artist(res_circle)

		t_cmap = "YlOrBr_r"
		ax = fig.add_subplot(gs[0,1])
		#pl(ax, press * m_air / (rho * kB), cmap=t_cmap, vbound=None)
		pl(ax, press / rho, cmap=t_cmap, vbound=[0,1.5])
		#pl(ax, press, cmap=t_cmap, vbound=None)

	gs.tight_layout(fig)
	plt.savefig("img/img_%05d.png" % fnum)
	plt.close()

pool = multiprocessing.Pool(8)
pool.map(doit, range(0, nfile))
