from dmhd.timestep import rk4_step
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
    