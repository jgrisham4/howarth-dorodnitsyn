#!/usr/bin/env python

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from dpath import path_dis

fs = 14
#rc('font',**{'family':'serif','serif':['Palatino']})
#rc('font', family='Times New Roman')
rc('text', usetex=True)
#rc('legend',**{'fontsize':fs,'numpoints':1})
#rc('axes',**{'labelsize':fs})

# Importing mach number data
fsize = (3.5, 3.5/1.33)
M = [1, 2, 3, 4]
fig1 = plt.figure(figsize=fsize)
markers = ["-k", "--k", "-.k", ":k"]
i = 0
ydata = []
Tdata = []
for Mnum in M:
    d = np.genfromtxt("results.dat".format(Mnum), skip_header=1)
    plt.plot(d[:, 1], d[:, 2], markers[i], label=r"$M_e = {}$".format(Mnum), lw=1.5)
    ydata.append(d[:, 1])
    Tdata.append(d[:, 3])
    i += 1
plt.xlabel(r"$\eta \, = \, \displaystyle y \sqrt{\frac{U_e}{\nu_e \, x}}$")
plt.ylabel(r"$u/U_e$")
plt.ylim([0.0, 1.1])
plt.xlim([0.0, 14.0])
plt.legend(loc=4)
plt.tight_layout()

fig2 = plt.figure(figsize=fsize)
for i in range(0, len(ydata)):
    plt.plot(ydata[i], Tdata[i], markers[i], label=r"$M_e = {}$".format(M[i]), lw=1.5)
plt.xlabel(r"$\eta \, = \, \displaystyle y \sqrt{\frac{U_e}{\nu_e \, x}}$")
plt.ylabel(r"$T/T_e$")
plt.ylim([1.0, 4.5])
plt.xlim([0.0, 14.0])
#plt.legend()
plt.tight_layout()

#fig1.savefig("../report/images/uqUe.pdf")
#fig2.savefig("../report/images/TqTe.pdf")

plt.show()
