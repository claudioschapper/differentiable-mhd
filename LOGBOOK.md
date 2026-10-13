# Logbook

## 2026-10-10 / 2026-10-11
Goal: 1D heat equation, pseudo-spectral, vs analytic solution.
- Spectral first derivative (`diffspec`) done: matches analytic to ~1e-14, Nyquist mode zeroed.
- Bugs hit: `N/2` float used as index, missing factor of i, negative wavenumbers had wrong sign, `if`/`if`/`else` instead of `elif`.
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
NEXT: stability experiment on the heat equation (h_crit = 2.785/(nu (N/2)^2); seed 1e-10 Nyquist noise, measure per-step growth vs R(-nu k^2 h)); then integrating-factor RK4 and 1D Burgers with/without 2/3 dealiasing.
