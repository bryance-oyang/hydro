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
rc("axes", grid=True)
rc("grid", linestyle="--")
rc("xtick", direction="in")
rc("ytick", direction="in")
rc("savefig", format="pdf", bbox="tight")
import numpy as np
import os

for fnum in range(200):
	print("fnum = %d" % fnum)

	if fnum == 0:
		param = np.loadtxt("data/param.dat")
		nx = param[0]
		ny = param[1]
		dx = param[2]
		dy = param[3]
	rho = np.transpose(np.loadtxt("data/rho_%05d.dat" % fnum))
	press = np.transpose(np.loadtxt("data/press_%05d.dat" % fnum))

	[x, y] = np.mgrid[0:(nx+1)*dx:(nx+1)*1j, 0:(ny+1)*dy:(ny+1)*1j]

	def pl(ax, q, cmap="viridis", vbound=None):
		if vbound is None:
			vmin = 0
			vmax = np.amax(q)
		else:
			[vmin, vmax] = vbound

		im = ax.pcolormesh(x.T, y.T, q, cmap=cmap, norm=colors.Normalize(vmin=vmin, vmax=vmax))
		ax.set_aspect(1)
		divider = make_axes_locatable(ax)
		cax = divider.append_axes("right", size="5%", pad=0.05)
		cbar = plt.colorbar(im, cax=cax)
		cbar.ax.tick_params()

	fig = plt.figure(figsize=(12, 6))
	gs = gridspec.GridSpec(1, 2)

	ax = fig.add_subplot(gs[0,0])
	pl(ax, rho, cmap="gray", vbound=[0,2.5])

	ax = fig.add_subplot(gs[0,1])
	pl(ax, press, cmap="gray", vbound=[0,3])

	gs.tight_layout(fig)
	#plt.show()
	plt.savefig("img/img_%05d.png" % fnum)
	plt.close()
