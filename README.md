# LIFEguidance

MATLAB implementation of the *weakly-constrained, fixed-phase, single-spacecraft
guidance problem at L2* of `Guidance_Formulation_LIFE.pdf` (Sec. 4). One LIFE
collector with 20 fixed, throttleable thrusters must reach a point fixed on a
rotating virtual formation near Sun–Earth L2 within t_f = 4 h, with minimum
total impulse and a free arrival time τ. The code structure follows the
earlier Malladi et al. (2019) reproduction on branch `malladi-et-al-claude`.

## Usage

Requires MATLAB (R2026a tested; Optimization and Symbolic toolboxes are needed
for the tests only) and an internet connection on first run. CasADi 3.8.1 with
IPOPT is downloaded into `external/` (git-ignored) by `setupCasadi.m`.

```matlab
main                 % solve + refine (~2 min), validate, figures, two videos (~10 min)
runtests('tests')    % model consistency checks (~5 s)
```

`main.m` has stage flags (`runSolve`, `runPlots`, `runVideos`). With
`runSolve = false`, it reloads `results/solution.mat`. All outputs are written
to `results/` (git-ignored).

## Method

- **Model** (`src/model/`)
  - Relative attitude: eqs. (10) and (13), with a Hamilton scalar-first quaternion for R^f_i.
  - Relative translation: eqs. (19)–(22), on the CR3BP of eqs. (14)–(15), with cannonball SRP on the spacecraft only.
  - Thruster wrench: eqs. (29)–(31).
  - The gravity difference uses Battin's f(q) form, so no cancellation occurs.
  - `relativeDynamics.m` is type-generic: the same code is used for the NLP (CasADi) and for the propagation (double).
- **Transcription** (`src/ocp/`)
  - Radau collocation with 3 points per interval and piecewise-constant thrust, solved with IPOPT/MUMPS.
  - The free arrival time enters via s = t/τ. On nonuniform meshes the objective is Σ h_n Σ_k T_k,n.
  - Scaling, initial guess and tolerances are as specified in Sec. 4.3.
  - In the NLP, the reference r_bf(t) is a degree-14 Chebyshev fit of the CR3BP halo on [0, t_f] (fit error 6e-5 m). It enters only through the tidal and SRP terms.
- **Mesh refinement** (`runMeshRefinement.m`)
  - The N = 200 uniform mesh is doubled to 400 (the spec's "doubling N" check). Refinement then bisects intervals with thrust switches, partial throttle, or local discretization error above 1e-3 of the arrival tolerances.
  - Each level is warm-started from the previous one.
  - A level is *accepted* when it meets all of the following:
    - IPOPT returns `Solve_Succeeded`,
    - the independent propagation passes eqs. (59)–(60),
    - J changes by less than 1%.
  - The reported solution is the last accepted level.
- **Validation** (`src/sim/`, `src/post/`)
  - Independent ode89 propagation of the optimized controls.
    - The relative state is co-integrated with the reference ODE (14), not the fit.
    - RelTol 1e-11 and AbsTol 1e-12 apply to the scaled states, and the integrator restarts at every control switch.
  - The same controls are also propagated through a separately written absolute model (spacecraft and formation origin in SER, inertial attitude). Its gate is 1% of the arrival tolerances.
- **Conventions**
  - R_x, R_y, R_z are active rotations. This is the only reading under which eq. (52) gives ω^b_cb = +Ω e₃ and eq. (53) gives R^c_f e₃ = +e₂.
  - Eq. (19) needs the SER-frame velocity difference, which is reconstructed as R^b_f(ṙ^f_pi + ω^f_bf × r^f_fi).

## Results (nominal case, `results/summary.json`)

| Quantity | Value |
|---|---|
| Final mesh | level 6, N = 770 (nonuniform), `Solve_Succeeded` |
| Arrival time τ | 14400 s, upper bound t_f active (multiplier 0.50) |
| Total impulse J | 131.569 N s (Δm = 0.061 kg at I_sp = 220 s) |
| Max thrust | 9.9992 N (T_max = 10 N), bound violation 0 |
| Terminal ‖r‖ / ‖ṙ‖ | 1.7e-7 m / 8.5e-11 m/s (gates 1e-2 m, 1e-5 m/s) |
| Terminal angle (eq. 60) / ‖ω‖ | 9.9e-8 rad (atan2 form 1.2e-8) / 3.7e-13 s⁻¹ (gates 1e-4 rad, 1e-7 s⁻¹) |
| Doubling N 200 → 400 | ΔJ = 0.56 % |
| Propagation vs collocation | ≤ 1.7e-7 m, 8.5e-11 m/s, 1.2e-8 rad, 6.2e-12 s⁻¹ |
| Absolute-model cross-check | ≤ 2.5e-8 m, 1.9e-12 m/s, 6.7e-13 rad (pass) |
| max \|‖q‖−1\| | 1.4e-8 collocation, 2.3e-15 propagation |

Mesh history (J in N s; the arrival checks are on the propagated state):

| Level | N | J | ΔJ | (59)–(60) |
|---|---|---|---|---|
| 0 | 200 | 132.781 | – | **fail** (θ = 2.3e-4 rad) |
| 1 | 400 | 132.045 | 5.6e-3 | pass |
| 2–5 | 426–612 | 131.77 → 131.58 | 2e-3 → 2e-4 | pass |
| 6 | 770 | 131.569 | 1.0e-4 | pass (**reported**) |
| 7 | 1066 | 131.564 | 3.4e-5 | IPOPT only "Solved To Acceptable Level" → not accepted, stop |

The solution has 7 burns: 4 in the first 870 s (a large reorientation combined
with the departure burn), one at about 13090 s, and 2 at arrival.

### Known limitations

- **Near-impulsive burns.** Burns 1, 5 and 7 halve in duration at every bisection, and their peak thrust has not saturated. Their shape is therefore mesh-dependent, which is reported as "features converged: 0".
  - Between levels 5 and 6, the inertial-frame force integrals of burns 3–7 change by ≤ 0.4%.
  - Burns 1 and 2 change by 2.4% and 1.2%; they trade impulse with each other.
- **Local optimum only.**
  - A cold start on a uniform N = 400 mesh converges to a different optimum, J = 138.08 N s.
  - In an earlier run, an over-refined level (N = 1369, not mesh-converged) reached J = 130.36 N s.
  - Better local optima therefore probably exist; a multi-start study is the natural next step.
- **IPOPT setting.** `bound_relax_factor = 0` is required. With IPOPT's default, `honor_original_bounds` clips the returned thrusts by ~1e-8·T_max. That leaves an ω bias of about 1e-9 s⁻¹ per interval, which grows to about 1e-4 rad of terminal attitude error (measured).

## Outputs (`results/`)

- **Figures**
  - `tracking_errors.png`: θ, ω, r, ṙ, with tolerance bands and magnified arrival views.
  - `control_wrench.png`: body-frame force and torque.
  - Thrust: `thrust_histories.png` (per pod), `thrust_burn_zoom.png`, `thrust_heatmap.png`.
  - `impulse.png`: cumulative impulse and per-thruster impulse.
  - `propagation_discrepancy.png`, `mesh_convergence.png`.
  - Static trajectories: `trajectory_formation.png`, `trajectory_SER.png`.
- **Videos**
  - `video_SER.mp4` and `video_formation.mp4` share one event-aware frame schedule.
  - Coasts play at ×600. Each burn plays in slow motion for at least 3 s, with a pause at each burn start.
  - Every frame shows the time, the playback rate, a 20-thruster activity panel and a burn timeline.
  - Plumes point opposite to the force, scaled with T/T_max. The spacecraft is drawn ×3.
- **Data**
  - `solution.mat` (v7.3) holds all levels, the reference, the propagation and the parameters.
  - `summary.json`, `settings.json`.
  - CSV files: thrust history, propagated and collocation states, mesh convergence.
  - `ipopt_level*.log`.

## Layout

- `src/lifeParameters.m`: every constant of Sec. 4.3, the thruster geometry, solver settings and gates.
- `src/model/`: frames, dynamics (relative, inertial-path check, absolute), state mappings.
- `src/ref/`: halo propagation and Chebyshev fit.
- `src/ocp/`: collocation NLP, initial guess, warm start, mesh refinement.
- `src/sim/`, `src/post/`: propagation, local errors, burns, validation, saving.
- `src/plot/`: figures and videos.
- `tests/ModelTest.m`: constants, thruster positive span, frame conventions, halo period, gravity difference against vpa, eqs. (19)–(22) against the inertial path, CasADi against double, terminal equilibrium, quaternion against DCM kinematics, absolute against relative propagation.
