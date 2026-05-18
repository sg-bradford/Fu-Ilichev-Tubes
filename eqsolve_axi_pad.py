import numpy as np
from config import *

# We seek a solution by assuming the truncated Fourier series for h(xi)

# We set a grid in xi
# NB: we exclude final point in period
# from collacation set

def eqsolve_axi_pad(Mx,Lx,h0_hat,u0_hat,ufar,cw,icount,w0):

    hx = 2*Lx/(2*Mx+1)
    x = np.linspace(-Lx,Lx-hx,2*Mx+1)
    kx = (np.pi/Lx)*np.hstack([np.array(range(Mx+1)), np.array(range(-Mx,0))])

    # Solve via Newton's method

    eps = 10**-6    # Numerical differentiation
    tol = 10**-8    # Covergence tolerance

    w = w0

    Nw = len(w)

    k = 1
    corr = 99
    while corr > tol:
        f = my_ODE(x, hx, Lx, w, kx)
    return -1



def my_ODE(x, hx, Lx, w, kx):

    # note that x is denoted Z in Fu & Il'ichev

    h_hat = np.zeros(Mx+1, dtype = complex)
    h_hat[0] = w[0] # zero frequency
    h_hat[1:int(Mx/2)+1] = w[1:int(Mx/2)+1]
    h_fft = np.hstack([h_hat, np.flip(np.conj(h_hat[1:Mx+1]))])       # POTENTIAL PROBLEM HERE, FFT ASYMMETRIC
    
    st = int(Mx/2)+2

    u_hat = np.zeros(Mx+1, dtype = complex)
    u_hat[0] = w[st - 1]    # zero frequency is zero since function u is odd
    u_hat[1:int(Mx/2)+1] = 1j*w[st:st + int(Mx/2)]
    u_fft = np.hstack([u_hat, np.flip(np.conj(u_hat[1:Mx+1]))])

    cw = w[-1]
    ufar = w[-2]

    hh = np.real(np.fft.ifft(h_fft))
    hh_x = np.real(np.fft.ifft(1j*kx * h_fft))
    hh_xx = np.real(np.fft.ifft((1j*kx)**2 * h_fft))

    uu = ufar*x + np.real(np.fft.ifft(u_fft))
    uu_x = ufar + np.real(np.fft.ifft(1j*kx * u_fft))
    uu_xx = np.real(np.fft.ifft((1j*kx)**2 * u_fft))

    # Elastic strain energy [Ghent model of Fu & Il'ichev eq. (4), (5)]

    zd = lam2inf + uu_x

    lam1 = hh/aa
    lam2 = np.sqrt(hh_x**2 + zd**2)

    # sig1, sig2 from /fu-dispersion/elastic

    lam1_e2 = lam1**2
    lam1_e4 = lam1**4
    lam2_e2 = lam2**2
    lam2_e4 = lam2**4

    lam_exp24 = lam1_e2 * lam2_e4
    lam_exp42 = lam1_e4 * lam2_e2

    sig1 = mu * (lam_exp42 - 1) * Jm/(-1 - lam_exp42 + (lam2_e2 * -lam2_e2 + Jm + 3) * lam1_e2)
    sig2 = mu * (lam_exp24 - 1) * Jm/(-1 - lam_exp42 + (lam2_e2 * -lam2_e2 + Jm + 3) * lam1_e2)

    # compute values at infinity

    sig1inf = mu*(lam1inf**4 * lam2inf**2 - 1) * Jm/(-1 - lam1inf**4 * lam2inf**2 + lam2inf**2 * (lam2inf**2 + Jm + 3) * lam1inf**2)
    sig2inf = mu*(lam1inf**2 * lam2inf**4 - 1) * Jm/(-1 - lam1inf**4 * lam2inf**2 + lam2inf**2 * (lam2inf**2 + Jm + 3) * lam1inf**2)

    P0 = sig1inf/(lam1inf**2 * lam2inf * aa)

    uterm = aa * sig2 * zd / lam2**2
    utermfft = np.fft.fft(uterm)
    uterm_x = np.real(np.fft.ifft(1j*kx * utermfft))

    wterm = aa * sig2 * hh_x / lam2**2
    wtermfft = np.fft.fft(wterm)
    wterm_x = np.real(np.fft.ifft(1j*kx * wtermfft))

    # Solve Fu & Il'ichev [2010] (2.6) in a travelling frame for pressure
    
    press = Eh * (-wterm_x + sig1 / lam1 + rhow * aa * cw * hh_xx) / (hh * zd)

    FF = (press - P0) / rho   # cf. eq. (3.4) in our [Parau's] AFM paper

    Sx = (lam2inf + uu_x) * hh_x

    print(Sx)
    return -1

