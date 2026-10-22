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

## 2026-10-19
Goal: close out timestep.py with a stiff check of `if_step` against plain `rk4_step`.
- `test_stability_and_blowup`: Bernoulli y' = -lam y + mu y^2 with lam = 70, mu = 1, dt = 0.05 (n_steps = 20), so lam*dt = 3.5, past RK4's real-axis limit of ~2.8 (|R(-3.5)| ~ 2.7 per step). Asserts `if_step` is within 1e-3 relative error of the exact solution and that `rk4_step` is more than 10x off.
- Measured: `if_step` relative error ~5e-4 at t = 3 dt. My first tolerance of 1e-10 failed. It assumed roundoff-level error, but this is truncation error of the nonlinear stages at lam*dt = 3.5, where the O(dt^4) constant grows with (lam dt)^p. The stepper was fine and the tolerance was wrong, so tolerances should come from a printed measurement, not a guess.
- Comparing against the exact value at T = 1 is meaningless here: y(1) ~ e^-70 ~ 4e-31, so any stepper that returns ~0 passes. Compare relative error at early times, where y is still O(1e-2..1).
- RK4's failure is partly the quadratic term, not only linear instability: once the amplified mode pushes y past lam/mu = 70, mu y^2 takes over. Written as `not (err <= 10)` so a nan from overflow still counts as failure.
- Limits of the test: the loop runs 4 iterations and the `i > 2` guard means both assertions fire once, at t = 3 dt. It shows accuracy and RK4 failure at one point, not boundedness to T = 1. The 1e-3 tolerance is only ~2x above the measured error; it is a regression guard, do not tighten it. n_steps = 20 fixes lam*dt = lam*T/n_steps, so changing it moves the stiff regime.
- 23 tests pass (`uv run pytest`).
NEXT: `burgers.py`. Periodic [0, 2 pi), u_t = -u u_x + nu u_xx, nu ~ 0.1, N ~ 64. Exact Cole-Hopf reference with theta(x, 0) = 1 + a cos x, 0 < a < 1 (start a ~ 0.5); derive u = -2 nu theta_x / theta at general t by hand and check periodicity and the PDE before coding. Check the t = 0 round trip through fft/ifft, plot |u_hat(k)| of the exact IC to confirm it reaches roundoff before the 2/3 cutoff. Then dt-convergence (slope 4), spatial convergence on log-linear axes, and the RK4-vs-IF stability figure on Burgers.

## 2026-10-20 .. 2026-10-22
Goal: 1D viscous Burgers (`solvers/burgers.py`), verified against an exact solution, as the Phase 0 deliverable.

### Advective term
- Conservative form: u u_x = (u^2/2)_x, so N_hat(k) = -i k * FFT[u^2/2](k), 2/3-masked. The product is formed in physical space and the derivative applied afterwards in Fourier space, so `ik` never multiplies `u` itself. This fits the `if_step` contract (N takes a physical field, returns a Fourier array) with one forward FFT per stage; the non-conservative u*u_x would need three extra transforms to get u_x from u alone.
- Mask and derivative are both diagonal in k and commute; the mask goes on the transform of the product. `.real` is taken before squaring so a drifting imaginary part cannot feed back.
- k[0] = 0 exactly, so N[0] = 0 for any field: N is a flux divergence, hence the mean of u is conserved to roundoff. The non-conservative form only conserves it approximately, since aliasing breaks the discrete equivalence.

### Exact reference (Cole-Hopf)
- u = -2 nu theta_x / theta solves Burgers exactly when theta solves the heat equation, so theta(x, 0) = 1 + a cos x evolves as theta = 1 + a cos x e^{-nu t} (mode k = 1, so the decay is e^{-nu k^2 t} = e^{-nu t}), giving u = 2 nu a sin x e^{-nu t} / (1 + a cos x e^{-nu t}).
- 0 < a < 1 keeps theta > 0 (the heat equation only damps), so u never blows up; u is the derivative of a periodic log, so its mean is exactly zero. a is the steepness knob: as a -> 1 theta nearly vanishes at x = pi and a front forms there, with width set by nu. Small-a limit: u ~ 2 nu a sin x e^{-nu t}, satisfying u_t = nu u_xx with the nonlinear term O(a^2).
- theta_x is differentiated BY HAND, not spectrally. If the reference used `diffspec`, an error in the derivative code would appear on both sides of the comparison and cancel. Same principle as z = 1/y in the Bernoulli test. Only u_exact(x, 0) and u_exact(x, T) are needed, not all t.
- The solver stays general (takes an `initial` array); only the reference is specific to theta_0 = 1 + a cos x.

### Diagnostics before any dynamics (tests/test_burgers.py)
- `test_u_exact`: finite-difference residual u_t + u u_x - nu u_xx at one point (h = 1e-4 first derivatives, 1e-3 second). Measured 7.23e-10, tolerance 1e-8. Finite differences, not analytic derivatives, so the test does not reuse the derivation it is checking.
- `test_N`: u = sin x gives -(u^2/2)_x = -1/2 sin 2x. Checks sign, the factor 1/2, the wavenumber layout, and that mode 1 feeds mode 2. Imaginary part also asserted zero.
- `test_N_mean`: field with a nonzero mean, N[0] == 0 exactly (not approximately).
- `test_exact_spectrum`: |u_hat[0]| at roundoff and |u_hat[n//3]|/max = 3.6e-12. The spectrum decays like k r^k with r = (1-sqrt(1-a^2))/a ~ 0.27 for a = 0.5, so at n = 64 the tail is ~1e-11 (NOT roundoff) and n = 96 is needed for the dt study to have no spatial floor.
- `test_aliasing`: u = cos(21x) (the top kept mode at n = 64). The product has content at k = 42, which aliases to -22, outside the mask. Correct framing: the 2/3 rule works on band-limited input, so the test feeds a mode INSIDE the mask and checks its alias is removed, rather than feeding a mode above n/3.
- `test_dynamics`: n = 96, nu = 0.1, a = 0.5, T = 1. Error 2.6e-15 at 200 steps, 1.5e-15 at 400: at the roundoff floor and not time-limited.

### dt-convergence (scripts/burgers_dt_convergence.py, figures/burgers_dt_convergence.png)
- n = 96, nu = 1, a = 0.5, T = 1, n_steps = 4..64. Errors 2.98e-5 .. 9.64e-10; observed orders 3.44, 3.71, 3.85, 3.92; fitted slope on the last three points 3.886. Approaches 4 from below, as in the Bernoulli test, and the finest point is far from roundoff so no plateau contaminates the fit.
- nu = 0.1 is useless as an order test: the amplitude is 2 nu a, so the problem is nearly linear, the integrating factor handles the linear part exactly, and every step count lands on the roundoff floor. nu = 1 raises the amplitude to 1 and is the first real test of the nonlinear stages. Larger errors at nu = 1 are therefore expected, not a regression.
- The spectrum depends only on `a`, not nu, so n = 96 keeps the spatial error negligible in both cases.

### Tolerance mistakes worth not repeating
- Every threshold here was set from a printed measurement, after three failures from guessed ones: 1e-10 on the stiff Bernoulli check (truncation error, not roundoff); 1e-17 then 1e-14 on the exact spectrum; 1e-14 on the aliasing noise (cos(21x) carries ~1e-14 of phase noise at argument ~132, the unnormalised FFT sums 64 of them, and the derivative scales by |k| up to 21, giving ~1e-13).
- Two test-shape bugs caught by inspection: `np.all(u_mask[22:43]) == 0` collapses the slice to one bool before comparing, so it only asserted that at least one slot was zero (it did still fail with the mask removed, but a mask zeroing one slot of 21 would have passed); and the second dynamics run called `solve_burgers` with `n_steps` instead of `2*n_steps`, so the step-doubling check tested nothing.
- Committed: f824e72 (burgers + tests), 05971cc (convergence script + figure).
NEXT: spatial convergence. Fix dt small enough that time error is below 1e-12 (nu = 0.1 with ~200 steps), vary n = 8, 12, 16, .., and plot max error vs n on log-LINEAR axes. Expect a straight falling line of slope ln r ~ -1.3 per mode for a = 0.5, flattening at ~1e-14 around n ~ 60. Then the IF-vs-RK4 stability figure on Burgers, and only then Phase 1 (2D grid, Poisson inversion with k = 0 explicit, Jacobian on an analytic pair, energy/enstrophy conservation at nu = 0).

## 2026-10-22 (spatial convergence)
Goal: spatial convergence of the Burgers solver, checked against a predicted decay rate rather than just "it goes down".
- Setup (scripts/burgers_spatial_convergence.py, figures/burgers_spatial_convergence.png): nu = 0.1, a = 0.5, T = 1, n_steps = 200 fixed so time error stays ~1e-15; n = 8, 16, .., 96. Plotted on log-LINEAR axes: spectral error is exponential in n, not a power of n, so there is no "order", only a decay rate.
- Prediction: if u is analytic in the strip |Im x| < y0, its Fourier coefficients decay like exp(-y0 |k|), and the truncation error at cutoff k_max scales like exp(-y0 k_max). u = -2 nu (ln theta)_x is singular where theta vanishes; with theta = 1 + a' cos x and x = pi + i y, cos(pi + i y) = -cosh y, so the zero is at cosh(y0) = 1/a' with a' = a e^{-nu t}. y0 is also the width of the front in x: a -> 1 sharpens it, y0 -> 0 and the required n diverges.
- Result: fitted decay rate 0.478 per unit n (natural log, fit on n >= 16 with errors above 1e-12). Prediction y0(T)/3 = 0.477. Ruled out: y0(0)/3 = 0.439 (8% off) and the full band y0/2 ~ 0.66-0.72.
- So the 2/3 mask sets the effective cutoff (k_max = n/3, not n/2): the solver resolves exactly the modes the dealiasing rule says it should.
- Correction to my own advice: I said to use y0 at t = 0, as the least smooth moment. Wrong. The error is measured at T, so it is the truncation tail of u(T). Two reasons: viscosity smooths the solution (a' shrinks, the singularity retreats from the real axis, y0 grows from 1.32 to 1.43); and truncation error made early sits at k ~ n/3 ~ 20, where the integrating factor damps it by e^{-nu k^2 t} = e^{-40 t}, so only recent error survives. The error is slaved to the current solution, not the initial condition.
- The third-digit agreement (0.2%) is partly luck: the true tail is ~ k r^k with a polynomial prefactor and the fit window shifts the number. The honest claim is agreement within ~5%, which is the test tolerance. It still rejects the t = 0 value and the full-band alternative.
- Floor at ~2.6e-15 for n >= 72, the same value as test_dynamics: roundoff, no information about the rate. Error at n = 64 is 1.3e-14.
- `test_spatial` (tests/test_burgers.py): asserts error at n = 64 < 1e-13, at least 4 points in the fit window (so polyfit cannot silently fit two), and |y0(T)/3 - rate|/rate < 5e-2. The n = 64 check sits after the loop so it fails loudly if 64 drops out of the grid.
- Committed: 4373725.
NEXT: stability figure (scripts/burgers_stability.py). Same Burgers problem, nu = 0.1, n = 64, so nu k_max^2 = 102.4 with k_max = n/2 (L acts on all modes, only N is masked). Choose dt with nu k_max^2 dt ~ 5 (about 20 steps to T = 1), past RK4's real-axis limit of 2.785. Plain RK4 on the full RHS L u_hat + N via rk4_step should blow up from a roundoff seed at |R(-5)| ~ 13.7 per step; IF-RK4 should stay accurate. Plot max|u_hat| vs t with the predicted growth line overlaid. Then Phase 1.

## 2026-10-22 (stability figure)
Goal: the figure that justifies the integrating factor — plain RK4 vs IF-RK4 on the same Burgers problem at the same step. Phase 0 deliverable complete.
- Setup (scripts/burgers_stability.py): nu = 0.1, a = 0.5, n = 64, T = 1, n_steps = 20 so dt = 0.05. k_max = 32 (L acts on all modes; only N is masked), so nu k_max^2 dt = 5.12, past RK4's real-axis stability limit of 2.785. Plain RK4 steps the full RHS L u_hat + N(ifft(u_hat)) through `rk4_step`; both runs use the SAME masked `N`, so the only difference between them is the treatment of L.
- figures/burgers_stability.png: amplitude and max-error vs t. IF-RK4 error flat at ~3e-11 for the whole run; plain RK4 climbs from roundoff and crosses the solution's size near t = 0.75, reaching 2.8e3 by t = 0.9.
- figures/burgers_stability_modes.png is the stronger result. Tracking |u_hat_k| per step for k = 32, 31, 30 gives measured growth 15.251, 11.460, 8.523 against predicted |R(-nu k^2 dt)| = 15.251, 11.460, 8.523 — three decimals, over 14 decades, with no bend. So the blow-up is EXACTLY the linear amplification of the grid-scale diffusion modes, not the nonlinearity and not aliasing.
- Why those lines never bend: the 2/3 mask kills N for |k| > 21, so modes 30-32 receive no nonlinear forcing and evolve under RK4 on the linear term alone at every amplitude. The solution blows up only because their huge u^2 content feeds back into the unmasked modes through `u.real`.
- The max-norm error fit (12.9-13.7 depending on the window) is NOT the number to quote: it mixes modes with different growth rates plus a pre-asymptotic plateau at RK4's own truncation error (~1e-10, same as IF-RK4's), so the fastest mode only dominates late. Per-mode tracking is the right diagnostic; the max-norm curve is the visual story.
- Mistake worth not repeating: my first blow-up guard tested max|u_hat| BEFORE the step and skipped stepping above 1e10. That cannot work here — the last step jumped 5.6e4 -> 3.8e35 in one step, because once u is large RK4 iterates the quadratic through its stages (k1 ~ u^2, k2 ~ u^4, k3 ~ u^8, k4 ~ u^16), so growth stops being ~15x per step and becomes a power law. A guard must test the RESULT of a step and discard it, not the input. The class of bug: a threshold on the input of an operation that amplifies super-exponentially.
- Two latent bugs found by inspection: the IF recording sat below a `continue`, so after RK4 blew up `errs_if` kept its np.zeros value and would have dropped to the log floor; and `ts = dt*np.arange(n_steps)` stopped at 0.95, so the run never reached T.
- The point, in scaling terms: plain RK4 needs dt < 2.785/(nu k_max^2) = 0.027 here (>= 37 steps), and that bound scales like 1/n^2 — at 256^2 in 2D the diffusion term alone would dictate the step. IF-RK4's limit is advective, dt ~ 1/(|u|_max k_max) ~ 1/n. That gap is why the integrating factor is in the solver.
- Committed: 73eb38b.
NEXT: Phase 0 is done (spectral derivative, RK4, heat, IF-RK4, Burgers with dt- and spatial-convergence and the stability figure). Phase 1: 2D grid and wavenumbers; Poisson inversion psi_hat = omega_hat/k^2 with the k = 0 mode set explicitly to zero; the Jacobian J(psi, omega) tested against an analytic pair; then nu = 0 with dealiasing must conserve energy and enstrophy to roundoff. Taylor-Green, then decaying turbulence with the k^-3 enstrophy cascade. Optional small test: mirror test_stability_and_blowup on the PDE (same dt, assert IF accurate and plain RK4 not, nan-safe).
