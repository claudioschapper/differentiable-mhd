import numpy as np
from dmhd.solvers.nse2d import solve_nse, nonlinear
from dmhd.timestep import if_step
from dmhd.spectral import Grid, partial_x, partial_y, laplacian, poisson_inversion, jacobian
import pytest


@pytest.mark.parametrize("N", [20, 24, 28])
@pytest.mark.parametrize("n_steps", [1, 10**2, 2*10**2, 10**3])
@pytest.mark.parametrize("l", [2*np.pi, 3])
def test_linear(monkeypatch, N, n_steps, l):
    monkeypatch.setattr(
        "dmhd.solvers.nse2d.nonlinear",
        lambda omega, l, grid: np.zeros((N, N), dtype=complex),
    )
    g = Grid.make(l,N)
    x, y = g.coords
    omega0 = np.sin(2*(2*np.pi/l)*x)*np.cos(3*(2*np.pi/l)*y)
    nu = 0.05
    T = 1
    omega_T = solve_nse(omega0, nu, T, n_steps, l)
    rate = 13 * nu * (2*np.pi/l)**2
    decay = np.exp(-rate * T)
    assert decay > 0.01
    assert np.max(np.abs(omega_T-omega0*decay)) < 1e-12

@pytest.mark.parametrize("N", [20, 24, 28])
@pytest.mark.parametrize("l", [2*np.pi, 3])
def test_nonlinear(N, l):
    g = Grid.make(l,N)
    x, y = g.coords

    omega = np.cos(3*(2*np.pi/l)*x+ 2*np.pi/l*y)+ np.cos(2*(2*np.pi/l)*x+ 2*np.pi/l*y)
    j_hat = nonlinear(omega, l, g)
    j_hat_hand = np.fft.fft2(-1/20*(np.cos((2*np.pi/l)*x) -np.cos(5*(2*np.pi/l)*x+2*(2*np.pi/l)*y)))
    max_j = np.max(np.abs(j_hat_hand))
    assert max_j > 0.1
    assert np.max(np.abs(j_hat-j_hat_hand)/max_j) < 1e-12
    
    predicted = np.zeros(j_hat.shape, dtype=bool)
    for i in [-1, 1]:
        predicted[i*5, i*2] = True
        predicted[i*1, 0] = True
    
    assert np.min(np.abs(j_hat[predicted])) > N**2/80
    assert np.max(np.abs(j_hat[~predicted])) / max_j < 1e-12



    