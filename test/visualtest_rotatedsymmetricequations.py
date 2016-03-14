import os
import numpy
import matplotlib.pyplot as plt

# plt.style.use("./paper-meanfieldmethod.mplstyle")

figsize = numpy.array([2, 1.5])*5.5
axesgeometry = [0.16, 0.18, 0.82, 0.78]

sourcedir = "."
# targetdir = "../images"


name = "rotsym.dat"
path = os.path.join(sourcedir, name)
x = numpy.loadtxt(path, delimiter=";")

T = x[:, 0]
sx_ind = x[:, 1]
sy_ind = x[:, 2]
sz_ind = x[:, 3]
sx_mf = x[:, 4]
sy_mf = x[:, 5]
sz_mf = x[:, 6]
sx_mfsym = x[:, 7]
sy_mfsym = x[:, 8]
sz_mfsym = x[:, 9]
sx_mpc = x[:, 10]
sy_mpc = x[:, 11]
sz_mpc = x[:, 12]
sx_master = x[:, 13]
sy_master = x[:, 14]
sz_master = x[:, 15]
error_ind = x[:, 16]
error_mf = x[:, 17]
error_mpc = x[:, 18]

name = "rotsym2.dat"
path = os.path.join(sourcedir, name)
x = numpy.loadtxt(path, delimiter=";")

T2 = x[:, 0]
sx_ind2 = x[:, 1]
sy_ind2 = x[:, 2]
sz_ind2 = x[:, 3]
sx_mf2 = x[:, 4]
sy_mf2 = x[:, 5]
sz_mf2 = x[:, 6]
sx_mfsym2 = x[:, 7]
sy_mfsym2 = x[:, 8]
sz_mfsym2 = x[:, 9]
sx_mpc2 = x[:, 10]
sy_mpc2 = x[:, 11]
sz_mpc2 = x[:, 12]
sx_master2 = x[:, 13]
sy_master2 = x[:, 14]
sz_master2 = x[:, 15]
error_ind2 = x[:, 16]
error_mf2 = x[:, 17]
error_mpc2 = x[:, 18]


fig_sx = plt.figure(figsize=figsize)
ax_sx = fig_sx.add_axes(axesgeometry)
#ax_sx.plot(T, sx_ind)
ax_sx.plot(T, sx_mf)
ax_sx.plot(T, sy_mf2)
ax_sx.plot(T, sx_mfsym, "x")
#ax_sx.plot(T, sx_mpc)
ax_sx.plot(T, sx_master, "--", color="b")
ax_sx.plot(T, sy_master2, "--", color="k")

ax_sx.set_xlim([0, 5])
ax_sx.set_xlabel(r"$\mathrm{Time}\ [\gamma^{-1}]$", labelpad=0.7)
ax_sx.set_ylim([-1, 1])
ax_sx.set_ylabel(r"$\langle\sigma_x\rangle$", labelpad=-5.0)


fig_sy = plt.figure(figsize=figsize)
ax_sy = fig_sy.add_axes(axesgeometry)
#ax_sy.plot(T, sy_ind)
ax_sy.plot(T, sy_mf)
ax_sy.plot(T, -sx_mf2)
ax_sy.plot(T, sy_mfsym, "x")
#ax_sy.plot(T, sy_mpc)
ax_sy.plot(T, sy_master, "--", color="b")
ax_sy.plot(T, -sx_master2, "--", color="k")

ax_sy.set_xlim([0, 5])
ax_sy.set_ylim([-1, 1])
ax_sy.set_xlabel(r"$\mathrm{Time}\ [\gamma^{-1}]$", labelpad=0.7)
ax_sy.set_ylabel(r"$\langle\sigma_y\rangle$", labelpad=-5.0)


fig_sz = plt.figure(figsize=figsize)
ax_sz = fig_sz.add_axes(axesgeometry)
#ax_sz.plot(T, sz_ind)
ax_sz.plot(T, sz_mf)
ax_sz.plot(T, sz_mf2)
ax_sz.plot(T, sz_mfsym, "x")
#ax_sz.plot(T, sz_mpc)
ax_sz.plot(T, sz_master, "--", color="b")
ax_sz.plot(T, sz_master2, "--", color="k")

ax_sz.set_xlim([0, 5])
ax_sz.set_ylim([-1, 1])
ax_sz.set_xlabel(r"$\mathrm{Time}\ [\gamma^{-1}]$", labelpad=0.7)
ax_sz.set_ylabel(r"$\langle\sigma_z\rangle$", labelpad=0.7)


# fig_error = plt.figure(figsize=figsize)
# ax_error = fig_error.add_axes(axesgeometry)
# ax_error.plot(T, error_ind)
# ax_error.plot(T, error_mf)
# ax_error.plot(T, error_mpc)

# ax_error.set_xlim([0, 5])
# ax_error.set_ylim([0, 1])
# ax_error.set_xlabel(r"$\mathrm{Time}\ [\gamma^{-1}]$", labelpad=0.7)
# ax_error.set_ylabel(r"Trace distance", labelpad=0.7)

plt.show()
