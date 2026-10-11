# Logbook

## 2026-10-10 / 2026-10-11
Goal: 1D heat equation, pseudo-spectral, vs analytic solution.
- Spectral first derivative (`diffspec`) done: matches analytic to ~1e-14, Nyquist mode zeroed.
- Bugs hit: `N/2` float used as index, missing factor of i, negative wavenumbers had wrong sign, `if`/`if`/`else` instead of `elif`.
- Tests: exactness over N and k, Nyquist, convergence on exp(sin x). All green.
- Convergence plot of exp(sin x) saved to figures/convergence_exp_sin.png.
NEXT: second derivative multiplier (-k^2, decide Nyquist), RK4 in timestep.py, heat equation vs exp(-nu k^2 t) sin(kx).
