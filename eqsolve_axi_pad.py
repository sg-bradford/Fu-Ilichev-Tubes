from elasticTube import *

# We seek a solution by assuming the truncated Fourier series for h(xi)

# We set a grid in xi
# NB: we exclude final point in period
# from collacation set

def eqsolve_axi_pad(Mx,Lx,h0_hat,u0_hat,bern,ufar,cw,icount,iload,w0):

    hx = 2*Lx/(2*Mx+1)
    x = np.linspace(-Lx,Lx-hx,2*Mx+1).T
    kx = (pi/Lx)*np.concatenate((np.array(range(Mx+1)), np.array(range(-Mx,0))))

    # Sort FT

    if iload != 1:
        h00 = h0_hat[0]