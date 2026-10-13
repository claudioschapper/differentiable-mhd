from dmhd.spectral import wavenumbers
from dmhd.timestep import rk4_step
import numpy as np
def solve_heat(initial, nu, T, n_steps):
    """initial is the initial state of the grid, while u_hat is the fourier transform/ fourier state representation thereof"""
    N = len(initial)
    h_t = T/n_steps
    wavenums = wavenumbers(N)
    u_hat = np.fft.fft(initial)
    def rhs(t, u_h):
        return -nu*(wavenums**2)*u_h
    for i in range(n_steps):
        u_hat = rk4_step(rhs, i*h_t, u_hat, h_t)
    u = np.fft.ifft(u_hat).real
    return u


    
    
    