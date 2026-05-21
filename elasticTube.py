import numpy as np
from config import *
import eqsolve_axi_pad
import matplotlib.pyplot as plt

pi = np.pi

# Computation of axisymmetric periodic waves on a ferrofluid column
# using the Fokas method

Nhamp = 1
hamp1 = 0.1
hamp2 = 0.1
if Nhamp == 1:
    dhamp = 0
else:
    dhamp = (hamp2-hamp1)/(Nhamp-1)

Nlam = 1
rlam1 = 13.0    # multiple of 2*pi
rlam2 = 13.0

if Nlam == 1:
    dlam = 0
else:
    dlam = 2*pi*(rlam2 - rlam1)/(Nlam - 1)

icount = 0

cww = np.array([Nhamp])
amp = np.array([Nhamp])
anorm = np.array([Nhamp])
uff = np.array([Nlam])
u1 = np.array([Nlam])
lamm = np.array([Nlam])

for ilam in range(Nlam):   # only executed once
    lam = rlam1*2*pi+ilam*dlam
    Lx = 0.5*lam
    kwave = 2*pi/lam

    # creating grid
    dxi = lam/(2*Mx + 1)
    xi = np.linspace(-Lx,Lx - dxi,2*Mx + 1)

    for iham in range(Nhamp):   # only executed once
        icount += 1
        ampl = hamp1 + iham*dhamp

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

        h0_hat = np.zeros(Mx + 1, dtype = complex)
        h0_hat[0] = wconv[0]
        h0_hat[1:int(Mx/2) + 1] = wconv[1:int(Mx/2) + 1]
        h0_fft = np.hstack([h0_hat, np.flip(np.conj(h0_hat[1:Mx + 1]))])
        hsol = np.real(np.fft.ifft(h0_fft))

        u0_hat = np.zeros(Mx + 1, dtype=complex)
        u0_hat[0] = 1j*wconv[st - 1]
        u0_hat[1:int(Mx/2) + 1] = 1j*wconv[st:st + int(Mx/2)]
        u0_fft = np.hstack([u0_hat, np.flip(np.conj(u0_hat[1:Mx + 1]))])
        ufar = wconv[-2]
        usol = ufar*xi + np.real(np.fft.ifft(u0_fft))
        kx = (np.pi/Lx)*np.hstack([np.array(range(Mx + 1)), np.array(range(-Mx, 0))])
        usold = ufar + np.real(np.fft.ifft(1j*kx*u0_fft))

        cww[iham] = wconv[-1]
        amp[iham] = ampl
        anorm[iham] = np.linalg.norm(hsol - aa, ord = np.inf)
        uff[ilam] = ufar
        u1[ilam] = usol[0]
        lamm[ilam] = lam

# some extra stuff

xp = np.linspace(-Lx, Lx, num = 100, endpoint = True)
fp = wconv[0]/(2*Mx + 1)
mfp = 1
for ip in range(1, int(Mx/2) + 1):
    jp = int(Mx/2) - ip + 1
    aco = (wconv[ip]*np.cos(ip*np.pi))/(2*Mx + 1)
    fp += aco*np.exp(1j*ip*2*np.pi*xp/lam) + np.conj(aco)*np.exp(-1j*ip*2*np.pi*xp/lam)

# Plotting solution

fig, ax = plt.subplots(3, 1)

ax[0].plot(xp, fp/mfp, '-r')
ax[0].plot(xi, hsol, 'ob')

ax[1].plot(xi, usol, '-k')
ax[1].plot(xi, usold, '-r')

ax[2].plot(xi, hsol, '-k')
ax[2].plot(xi, 1 - hsol, '-r')

plt.show()






