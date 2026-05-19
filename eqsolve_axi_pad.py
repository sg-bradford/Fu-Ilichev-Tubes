import numpy as np
import scipy as sp
from config import *

# We seek a solution by assuming the truncated Fourier series for h(xi)

# We set a grid in xi
# NB: we exclude final point in period
# from collacation set

def eqsolve_axi_pad(Mx,Lx,h0_hat,u0_hat,ufar,cw,icount,w0, ampl):

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
        f = my_ODE(x, hx, Lx, w, kx, ampl)
    
        # build Jacobian

        jac = np.zeros((Nw, Nw))

        rat = 1
        for j in range(Nw):
            w[j] += eps * rat
            f1 = my_ODE(x, hx, Lx, w, kx, ampl)
            w[j] -= eps * rat

            for i in range(Nw):
                jac[i, j] = (f1[i] - f[i]) / eps

        # solve linear Jacobian system

        rhs = -f
        adj = np.linalg.solve(jac, rhs)
        djac = abs(np.linalg.det(jac))

        # update guesses

        wold = w
        ssum = 0
        for i in range(Nw):
            ssum += rhs[i] ** 2
            w[i] += adj[i]

        corr = np.sqrt(abs(ssum))
        print("     k       corr")
        print(str(k) + "    " + str(corr))

        if k > 50:
            print("******** Newton iterations did not converge! *********")
        
    w = wold
    print("Converged!")




def my_ODE(x, hx, Lx, w, kx, ampl):

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

    sig1 = mu * (lam_exp42 - 1) * Jm/(-1 - lam_exp42 + lam2_e2 * (-lam2_e2 + Jm + 3) * lam1_e2)
    sig2 = mu * (lam_exp24 - 1) * Jm/(-1 - lam_exp42 + lam2_e2 * (-lam2_e2 + Jm + 3) * lam1_e2)

    # compute values at infinity

    sig1inf = mu*(lam1inf**4 * lam2inf**2 - 1) * Jm/(-1 - lam1inf**4 * lam2inf**2 + lam2inf**2 * (lam2inf**2 + Jm + 3) * lam1inf**2)
    # unused variable - could be a problem?
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

    pp = np.sqrt((1 + Sx ** 2) * ((lam2inf * cw - uinf) ** 2 - 2 * FF))

    npp = np.linalg.norm(pp)

    # Integrand of Fokas Bernoulli integral

    fokint = np.zeros(int(Mx / 2))
    eulerian_period = 2 * Lx * (lam2inf + ufar)

    for nn in range(1, int(Mx / 2) + 1):
        
        kk = 2 * np.pi * nn / eulerian_period
        ii1 = sp.special.iv(1, kk * hh)
        mk = max(hh)
        sc = sp.special.iv(1, kk * mk)
        xp = lam2inf * x + uu   # Convert Eulerian to Lagrangian (see p. 4 of Notes)
        dxp = lam2inf + uu_x
        fokas = (1 / sc) * hh * pp * ii1 * np.cos(kk * xp) * dxp    # assume wave even in x

        # Evaluate Fokas integral using periodic trapezium rule

        fokint[nn - 1] = hx * sum(fokas)

    # Fu & Il'ichev [2010] equation (2.5) in a travelling frame
    fu = uterm_x - (press / Eh) * hh * hh_x - rhow * aa * cw ** 2 * uu_xx
    fuill = np.imag(np.fft.fft(fu))

    # Bernoulli condition (do *not* include k = 0)

    F_col = np.zeros(Mx + 4)

    F_col[:int(Mx / 2)] = fokint[:int(Mx / 2)]
    F_col[int(Mx / 2):Mx] = fuill[1:int(Mx / 2) + 1]
    F_col[Mx] = w[int(Mx / 2) + 1]
    
    F_col[Mx + 1] = hh[0] - lam1inf * aa
    F_col[Mx + 2] = (max(hh) - min(hh)) - 2 * ampl
    F_col[Mx + 3] = uu_x[0]

    return F_col

