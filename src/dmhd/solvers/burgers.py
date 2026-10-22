from dmhd.spectral import wavenumbers
from dmhd.timestep import if_step
import numpy as np

def solve_burgers(initial, nu, T, n_steps):
    """Evolve u_t + u u_x = nu u_xx on [0, 2 pi) (periodic) to time T.

    initial: real physical-space field on the grid (length n).
    Returns the real physical-space field at t = T, using n_steps of if_step.
    """
    n = len(initial)
    h_t = T/n_steps
    wavenums = wavenumbers(n)
    u_hat = np.fft.fft(initial) 
    L = -nu*(wavenums**2)
    for i in range(n_steps):
        u_hat = if_step(u_hat, L, lambda u: N(n, u), h_t)
    u = np.fft.ifft(u_hat).real
    return u

def u_exact(x, t, nu, a):
    """Exact Burgers solution from Cole-Hopf with theta(x, 0) = 1 + a cos x (0 < a < 1).

    u = -2 nu theta_x / theta. Closed form, no time-stepping and no FFT,
    so it is an independent reference for the solver.
    """
    return 2*nu*a*np.sin(x)*np.exp(-nu*t)/(1+a*np.cos(x)*np.exp(-nu*t))

def N(n, u):
    """Nonlinear term -(u^2/2)_x in Fourier space, 2/3-dealiased.

    u: physical-space field of length n (the real part is taken before squaring).
    Returns a length-n Fourier array: -ik * FFT[u^2/2], masked to |k| <= n//3.
    Matches the if_step contract (N takes a physical field, returns a Fourier array).
    """

    wavenums = wavenumbers(n)
    k_masked = np.zeros(n)
    k_masked[:n//3+1], k_masked[(2*n+2)//3:] = wavenums[:n//3+1], wavenums[(2*n+2)//3:]
    return (-1j*k_masked*np.fft.fft((u.real)**2)/2)