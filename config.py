# global vars

Nmode = 64  # desired number of Fourier modes

Mx = 128    # We use nMode modes and pad the Fourier series of h(x)
# with nMode zeros. This accounts for aliasing errors which cannot
# distinguish between exp(i*N*t) and exp(i*2*N*t) in the trapezium rule.
# See Trefethen & Weideman SIAM Review 56 (2014), fig. 3.1.

aa = 1.0    # undisturbed radius of elastic tube [R in Fu & Ill'ichev]
rho = 1.0   # fluid density
rhow = 1.0  # membrane density
mu = 1.0    # stain energy parameter
Jm = 1.0    # another strain energy parametermatrices in numpy
Eh = 1.0    # elastic parameter
lam1inf = 1.0
lam2inf = 1.0
uinf = 0.0
