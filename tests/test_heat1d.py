import numpy as np
from dmhd.solvers.heat1d import solve_heat
import pytest

def test_error():
    T = 1
    v = 0.1
    N = 32
    n_steps = 1000
    h = 2*np.pi/N
    x = h*np.arange(N)
    initial = np.sin(x) + 0.5*np.sin(3*x)
    sol = solve_heat(initial, v, T, n_steps)
    analytic = np.exp(-v*T)*np.sin(x) + 0.5*np.exp(-9*v*T)*np.sin(3*x)
    assert np.max(np.abs(sol-analytic)) < 1e-12

@pytest.mark.parametrize("c", [1, -1])
def test_mean(c):
    T = 1
    v = 0.1
    N = 32
    n_steps = 1000
    h = 2*np.pi/N
    x = h*np.arange(N)
    initial = np.sin(x) + 0.5*np.sin(3*x) + c
    sol = solve_heat(initial, v, T, n_steps)
    assert np.abs(np.mean(sol)-c) < 1e-12

# @pytest.mark.parametrize("n_steps", [5, 10, 20, 40])
def test_order():
    steps = [5, 10, 20, 40]
    T = 1
    v = 0.1
    N = 32
    hs = np.array([T/i for i in steps])
    x = 2*np.pi/N*np.arange(N)
    initial = np.sin(x) + 0.5*np.sin(3*x)
    analytic = np.exp(-v*T)*np.sin(x) + 0.5*np.exp(-9*v*T)*np.sin(3*x)
    errs = np.zeros(len(steps))
    for a, n_steps in enumerate(steps):
        sol = solve_heat(initial, v, T, n_steps)
        errs[a] = np.max(np.abs(sol-analytic))
    slope, _ = np.polyfit(np.log(hs), np.log(errs), 1)
    assert abs(4-slope) < 0.15