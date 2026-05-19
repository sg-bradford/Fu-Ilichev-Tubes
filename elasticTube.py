import numpy as np
from config import *
import eqsolve_axi_pad

pi = np.pi

# Computation of axisymmetric periodic waves on a ferrofluid column
# using the Fokas method

Nhamp = 1
hamp1 = 0.1
hamp2 = 0.1
if Nhamp == 1:
    dhamp = 0
else:
    dhamp = (hamp2 - hamp1)/(Nhamp-1)

Nlam = 1
rlam1 = 13.0    # multiple of 2*pi
rlam2 = 13.0

if Nlam == 1:
    dlam = 0
else:
    dlam = 2*pi*(rlam2-rlam1)/(Nlam-1)

icount = 0

w0 = []

for ilam in range(1, Nlam+1):   # only executed once
    lam = rlam1*2*pi + (ilam-1)*dlam
    Lx = 0.5*lam
    kwave = 2*pi/lam

    # creating grid
    dxi = lam/(2*Mx+1)
    xi = np.linspace(-Lx,Lx-dxi,2*Mx+1)

    for iham in range(1, Nhamp+1):
        icount += 1
        ampl = hamp1 + (iham-1)*dhamp

    # making an initial guess

    if icount == 1:
        w0 = np.array(np.genfromtxt('wsave.csv', delimiter=''))
        Mxsave = int(np.genfromtxt('Mxsave.csv', delimiter=''))
        cwave = []
        hh0_hat = []
        uu0_hat = []
        ufar = []

    # Main Calculation

    wconv = eqsolve_axi_pad.eqsolve_axi_pad(Mx, Lx, hh0_hat, uu0_hat, ufar, cwave, icount, w0, ampl)

    h0_hat = np.zeros(Mx+1, dtype = complex)
    h0_hat[0] = wconv[0]

    



