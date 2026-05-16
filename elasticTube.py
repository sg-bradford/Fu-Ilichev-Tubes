import numpy as np
import csv

pi = np.pi

# Computation of axisymmetric periodic waves on a ferrofluid column
# using the Fokas method

aa = 1.0    # undisturbed radius of elastic tube [R in Fu & Ill'ichev]
rho = 1.0   # fluid density
rhow = 1.0  # membrane density
mu = 1.0    # stain energy parameter
Jm = 1.0    # another strain energy parametermatrices in numpy
Eh = 1.0    # elastic parameter
lam1inf = 1.0
lam2inf = 1.0
uinf = 0.0


l1i = lam1inf
l2i = lam2inf
R = aa

Nmode = 64  # desired number of Fourier modes

Mx = 2*Nmode    # We use nMode modes and pad the Fourier series of h(x)
# with nMode zeros. This accounts for aliasing errors which cannot
# distinguish between exp(i*N*t) and exp(i*2*N*t) in the trapezium rule.
# See Trefethen & Weideman SIAM Review 56 (2014), fig. 3.1.

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

wconv = []
w0 = []

for ilam in range(1, Nlam+1):
    lam = rlam1*2*pi + (ilam-1)*dlam
    Lx = 0.5*lam
    kwave = 2*pi/lam

    # creating grid
    dxi = lam/(2*Mx+1)
    xi = np.linspace(-Lx,Lx-dxi,2*Mx+1).reshape(-1, 1)

    for iham in range(1, Nhamp+1):
        icount += 1
        ampl = hamp1 + (iham-1)*dhamp

    # making an initial guess

    if icount == 1:
        bern = -99  # irrelevant in this code

        w0 = np.array(np.genfromtxt('wsave.csv', delimiter='')).reshape(-1, 1)
        Mxsave = int(np.genfromtxt('Mxsave.csv', delimiter=''))
        cwave = []
        hh0_hat = []
        uu0_hat = []
        ufar = []

        if Mx > Mxsave:
            cwsave = w0[-1]
            ufar = w0[-2]
            Mdiff = (Nmode-1) - (Mxsave/2+1)
            fr = (2*Nmode+1)/(Mxsave+1)
            w0 = [fr*w0[:-1],np.zeros((Mdiff,1)),cwsave]

        if Mx < Mxsave:
            cwsave = w0[-1]
            ufar = w0[-2]
            w0 = [w0[:Nmode+1],cwsave]

    if icount > 1:
        w0 = wconv

    # Main Calculation





