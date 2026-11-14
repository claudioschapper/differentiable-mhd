from dmhd.spectral import Grid, jacobian, poisson_inversion
from dmhd.timestep import if_step
import numpy as np

def solve_nse(omega0, nu, T, n_steps, l):
    """Evolve 2D incompressible Navier-Stokes in vorticity form on a periodic
    [0, l)^2 box to time T:

        d_t omega = J(psi, omega) + nu * laplacian(omega)

    Conventions:
        omega = d_x v_y - d_y v_x            (vorticity)
        omega = -laplacian(psi)              (Poisson relation)
        v = (d_y psi, -d_x psi)              (divergence-free velocity)
        J(psi, omega) = psi_x omega_y - psi_y omega_x = -(v . grad omega)

    Assumes zero mean flow (k = 0 mode of v is fixed to zero; it is conserved
    and removable by a Galilean shift on a periodic box). No forcing, no friction.

    omega0: real physical-space vorticity on an (n, n) grid.
    Returns the real physical-space vorticity at t = T after n_steps of if_step.
    """
    n = len(omega0[0])
    assert omega0.shape == (n,n)
    g = Grid.make(l, n)
    h_t = T/n_steps
    omega_hat = np.fft.fft2(omega0) 
    
    L = -nu*(g.ky_even_der**2 +g.kx_even_der**2)
    for i in range(n_steps):
        # RHS = L*omega_hat + N(omega): L = -nu|k|^2 (diffusion, integrating factor),
        # N = FFT[+J(psi, omega)] since v.grad(omega) = -J with v = (psi_y, -psi_x)
        omega_hat = if_step(omega_hat, L, lambda u: nonlinear(u, l, g), h_t)

    return np.fft.ifft2(omega_hat).real


def nonlinear(omega, l, grid):
    """Nonlinear term J(psi, omega), 2/3-dealiased. (coming from advection term)

    u: physical-space field on n x n grid (the real part is taken before squaring).
    Returns a Fourier array: FFT2[J(poission_inversion(omega), omega)], masked to |k| <= n//3.
    Where n is the number of grid points per axis, l is the length of the grid.
    Matches the if_step contract (N takes a physical field, returns a Fourier array).
    """
    return jacobian(poisson_inversion(omega, grid), omega, grid)