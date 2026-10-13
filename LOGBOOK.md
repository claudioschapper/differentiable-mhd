# Logbook

## 2026-10-10 / 2026-10-11
Goal: 1D heat equation, pseudo-spectral, vs analytic solution.
- Spectral first derivative (`diffspec`) done: matches analytic to ~1e-14, Nyquist mode zeroed.
- Bugs hit: `N/2` float used as index, missing factor of i, negative wavenumbers had wrong sign.
- Tests: exactness over N and k, Nyquist, convergence on exp(sin x). All green.
- Convergence plot of exp(sin x) saved to figures/convergence_exp_sin.png.
NEXT: second derivative multiplier (-k^2, decide Nyquist), RK4 in timestep.py, heat equation vs exp(-nu k^2 t) sin(kx).

## 2026-10-13
Goal: 1D heat equation, pseudo-spectral + RK4, verified against the analytic solution.
- `diffspec_gen(v, m)`: general order-m derivative. Nyquist slot zeroed for odd m only. For even m, (ik)^m is even in k, so +N/2 and -N/2 agree and nothing is biased; keeping the mode is correct. Test: (-1)^j under m=2 returns -(N/2)^2 v.
- `rk4_step(f, t, y, h)` in timestep.py. Verified on y' = lam*y (lam = -1..-5): fitted order 4.02-4.12. On y' = t^3 the error is at roundoff (<1e-16), so the stage times are right.
- `solve_heat` (evolves u_hat, diagonal RHS -nu k^2 u_hat): matches e^(-nu t) sin x + 0.5 e^(-9 nu t) sin 3x below 1e-12 (n_steps=1000, N=32). Mean conserved, imaginary part of the final ifft ~5e-17. Fitted time order 4.06.
- Bug hit: order fit gave slope -4 because I regressed against n_steps/N instead of h = T/n_steps. A slope of exactly -p means the x-axis is inverse to the step.
- Tooling lesson: editing files on disk while VS Code held unsaved buffers caused a save conflict and an overwritten LOGBOOK. Save/close tabs first.
- RK4 stability verified on heat (notebooks/2026-10-13-heat-stability.ipynb): measured per-step Nyquist growth = R(z) to 7e-16; h_crit = 2.7853/(nu (N/2)^2) = 0.1088 (N=32, nu=0.1), growth there = 1 to ~1e-14. Blow-up runs match eps*R^n with the same roundoff seed eps ~ 7e-17; h = 0.111 (just above h_crit) only reaches 5e-3 by T=40, so "bounded at finite T" does not prove stability.
- test_timestep now asserts y' = 10y matches e^10; mean-conservation tolerance loosened to 1e-12.
NEXT: integrating-factor RK4 in timestep.py (h_crit ~ 1/(nu N^2) is the stiffness it removes); then 1D Burgers with/without 2/3 dealiasing.
