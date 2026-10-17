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

## 2026-10-17
Goal: integrating-factor RK4 (`if_step` in timestep.py), verified before any Burgers work. Source: Trefethen, Spectral Methods in MATLAB, Ch. 10 (Program 27); Kassam & Trefethen 2005 is the ETDRK4 follow-up, not read yet.
- Derivation (interaction picture): for u_hat_t = L u_hat + N(u_hat), set U = e^{-Lt} u_hat. Then U_t = e^{-Lt} N(e^{Lt} U), so L drops out of the linear part and RK4 only sees N. Stages need e^{L dt/2} and e^{L dt}. Convention chosen: L on the right-hand side, E = exp(L dt/2) with entries <= 1 (no positive exponentials, so no overflow for diffusion). Physical meaning: E is the exact viscous damping over half a step, so dt is limited by the nonlinear timescale (advective CFL), not by 1/(nu k_max^2).
- Derivation errors caught on the way: missing dt factors, wrong sign in the exponentials, missing the factor 2 on the middle RK4 weights, reusing a stage value instead of the new one.
- Implementation checks (tests/test_timestep.py):
  - `NoNonlinear`: N = 0, nu k_max^2 dt ~ 1e3 (plain RK4 would blow up) -> matches exp(L dt) u_hat to 1e-14.
  - `NoLinear`: L = 0, N = w^2 -> agrees with `rk4_step` to 1e-14 over 100 steps (this caught a missing N() call on stage a4).
  - `test_order`: Bernoulli y' = -lam y + mu y^2, lam=5, mu=1, y0=1, T=1, exact y = 1/(mu/lam + (1 - mu/lam) e^{lam}) from z = 1/y (z' = lam z - mu, linear). This is the only test that exercises the E-powers on the nonlinear stages, since the other two have E = 1 or no stages. Error ratios per halving give observed order 3.83, 3.92, 3.96, 3.98 (n = 10..160), approaching 4 from below. The 5% tolerance on log-ratio passes at the coarsest pair with little margin (3.83 -> ratio 14.2); do not tighten it.
- Bugs hit: `np.max` on a complex array compares lexicographically, not by magnitude, so it silently ignored most of the error (use `np.abs` first); in the NoLinear loop `u` was overwritten before being compared with the `rk4_step` output (stale-state comparison).
- Not done yet: stiff check in the Bernoulli problem (lam=100, dt=0.05, lam*dt=5 beyond RK4's ~2.8 limit; assert bounded and accurate for `if_step`, unstable for `rk4_step`; no order assertion since order reduction at large lam*dt is expected). This is the miniature of the logbook figure.
- Committed: bf1980a (timestep: integrating-factor RK4 with linear, nonlinear and order tests).
NEXT: stiff Bernoulli check; `burgers.py` (L = -nu k^2 on the layout of `wavenumbers(n)`, dealias the product with the 2/3 mask, take .real before squaring, Cole-Hopf initial condition checked against the setup before running dynamics); dt-convergence (slope 4), then spatial convergence on log-linear axes. No 2D until that plot is exponential.
