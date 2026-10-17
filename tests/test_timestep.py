from dmhd.timestep import rk4_step
from dmhd.timestep import if_step
from dmhd.spectral import wavenumbers
import numpy as np
import pytest




def test_timestep():
    T = 1
    n_steps = 10**4
    h = T/n_steps
    y = 1
    const = 10
    analytic = np.exp(10*T)
    def rhs(t, y):
        return const*y
    
    for i in range(n_steps):
        y = rk4_step(rhs, 0, y, h)
    
    assert(np.abs(np.float64(y)-analytic)/np.abs(analytic)) < 1e-13


def test_IFNoNonlinear():
    rng = np.random.default_rng(0)  
    n = 32
    w = wavenumbers(n)
    v = 0.01
    L = -v*(w**2)
    # want v * (k_max)**2 * dt approx 1e3 to test stability of stepper on the linear part RK4 would be unstable due to high wavenumber stiffness
    N = lambda w: np.zeros(w.shape, complex)
    dt = 400
    u_hat = rng.standard_normal(n) + 1j * rng.standard_normal(n)
    assert np.max(np.abs(np.exp(L*dt)*u_hat - if_step(u_hat, L, N, dt))) < 1e-14

def test_IFNoLinear():
    rng = np.random.default_rng(0)  
    n_steps = 10**2
    T = 1
    dt = T/n_steps
    u = rng.standard_normal(1) + 1j * rng.standard_normal(1)
    v = u
    L = np.zeros(1)
    N = lambda w: w**2
    for i in range(n_steps):
        u = if_step(u, L, N, dt)
        v = rk4_step(lambda t, y: y**2, i*dt, v, dt)
        assert np.abs(v - u) < 1e-14

def test_order():
    """uses the bernoulli equation to test the powers of E in the nonlinear stages"""
    T = 1
    lam = 5
    mu = 1
    L = np.array([-lam])
    N = lambda w: mu*w**2
    analytic = np.array([1/(mu/lam+(1-mu/lam)*np.exp(lam))])
    n_steps = [10, 20, 40, 80]
    errs = np.zeros(len(n_steps))
    for a, steps in enumerate(n_steps):
        y = np.array([1], complex)
        dt = T/steps
        for i in range(steps):
            y = if_step(y, L, N, dt)
        err = np.abs(analytic-y)
        if a > 0:
            assert np.abs((np.log(preverr/err) - np.log(16))/np.log(16)) < 0.05
        preverr = err
    