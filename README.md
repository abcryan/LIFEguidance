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
main                 % solve, validate, save results/ and figures (~45 s)
makeVideo            % render the two videos from the saved run (~1 min)
runtests('tests')    % model consistency checks (~5 s)
```

## Method

- **Model** (`src/model/`)
  - Relative attitude: eqs. (10) and (13), with a Hamilton scalar-first quaternion for R^f_i.
  - Relative translation: eqs. (19)–(22), on the CR3BP of eqs. (14)–(15), with cannonball SRP on the spacecraft only.
  - Thruster wrench: eqs. (29)–(31).
  - The gravity difference uses Battin's f(q) form, so no cancellation occurs.
  - `relativeDynamics.m` is type-generic: the same code is used for the NLP (CasADi) and for the propagation (double).
- **Transcription** (`src/ocp/`)
  - One uniform mesh of `p.ocp.N` = 400 intervals.
  - Radau collocation with 3 points per interval and piecewise-constant thrust, solved with IPOPT/MUMPS.
  - The free arrival time enters via s = t/τ, and the objective is eq. (58).
  - Scaling, initial guess and tolerances are as specified in Sec. 4.3.
  - In the NLP, the reference r_bf(t) is a degree-14 Chebyshev fit of the CR3BP halo on [0, t_f] (fit error 6e-5 m). It enters only through the tidal and SRP terms.
- **Validation** (`src/sim/`, `src/post/`)
  - Independent ode89 propagation of the optimized controls.
    - The relative state is co-integrated with the reference ODE (14), not the fit.
    - RelTol 1e-11 and AbsTol 1e-12 apply to the scaled states, and the integrator restarts at every control switch.
  - The same controls are also propagated through a separately written absolute model (spacecraft and formation origin in SER, inertial attitude). Its gate is 1% of the arrival tolerances.
  - All reported errors, figures and videos are based on the independent propagation.
- **Conventions**
  - R_x, R_y, R_z are active rotations. This is the only reading under which eq. (52) gives ω^b_cb = +Ω e₃ and eq. (53) gives R^c_f e₃ = +e₂.
  - Eq. (19) needs the SER-frame velocity difference, which is reconstructed as R^b_f(ṙ^f_pi + ω^f_bf × r^f_fi).

## Virtual formations and initial condition

Both formations are built the same way, from a spec in `lifeParameters.m` (`p.target`, `p.init`):

- **Shape.** A pyramid with four collectors at the corners of a nulling (x) × imaging (y) rectangle in its z = 0 plane, centred on the frame origin, and the combiner at the apex, `height` along +z.
  - Points 1–4 are the corners (+,+), (−,+), (−,−), (+,−) in (x, y). Point 5 is the apex.
- **Attitude.** The frame points at the `sky` direction (a unit vector in the inertial frame) and spins about it: R^c = skyBasis(sky)·R_z(phase0 + rate·t).
  - skyBasis takes x = sky × e₃. With sky = e₂ this reproduces eq. (53).

| | Target formation | Initial formation (exists only at t = 0) |
|---|---|---|
| Rectangle (nulling × imaging) | 80 × 480 m | 25 × 150 m |
| Apex height | 960 m | 300 m |
| Sky direction | e₂ | −e₂ (opposite) |
| Spin | 2π/24 h about its own z | 2π/24 h about its own z (the same line, opposite axis) |
| Point used | 1 → target point r^f_fp = [40, 240, 0] m | 1 → spacecraft start |

The initial formation shares the target origin at t = 0, in both position and velocity. The spacecraft starts at the chosen point of the initial formation, with its body axes aligned with the initial frame, and co-rotates with it. That means:

- R^c_i(0) = R^c_f0(0);
- ω^i_ci(0) = rate·e₃;
- the inertial velocity relative to the origin is ω_cf0 × r.

This replaces the fixed initial state of eqs. (55)–(56). See `initialRelativeState.m`.

For the defaults, the initial attitude error is 180° about y_f, and the start lies 173 m from the target corner. To change the setup, edit the sky direction, baselines, apex height, spin rate or phase, or the point, then call `deriveFormations(p)`.

`tests/FormationTest.m` checks the construction. The central test maps the initial state into the absolute model and compares attitude, body rate, position and inertial velocity with the initial pyramid. It does this for all 5 points, for the default and for a randomly oriented, non-default initial formation.

## Results (nominal case, N = 400, T_max = 3 N, `results/summary.json`)

T_max is set to 3 N in `lifeParameters.m`; the document specifies 10 N.

| Quantity | Value |
|---|---|
| Solver | `Solve_Succeeded`, 212 iterations, 20 s (`main` total ~45 s) |
| Arrival time τ | 14400 s, upper bound t_f active (multiplier 0.72) |
| Total impulse J | 174.266 N s (Δm = 0.081 kg at I_sp = 220 s) |
| Max thrust | 1.66 N (T_max = 3 N, bound not active), bound violation 0 |
| Terminal ‖r‖ / ‖ṙ‖ | 1.3e-6 m / 1.4e-9 m/s (gates 1e-2 m, 1e-5 m/s) |
| Terminal angle (eq. 60) / ‖ω‖ | below the acos round-off floor (atan2 form 7.4e-8 rad) / 1.5e-11 s⁻¹ (gates 1e-4 rad, 1e-7 s⁻¹) |
| Absolute-model cross-check | ≤ 1.7e-8 m, 1.2e-12 m/s, 1.5e-13 rad (pass) |
| max \|‖q‖−1\| | 1.6e-6 collocation, 4.0e-15 propagation |

Choice of N (measured with the document's initial state and T_max = 10 N):
with N = 200, all checks pass except the attitude check (60), which gives
θ = 2.3e-4 rad.

### Known limitations

- **Burns are not resolved.** The optimum is near-impulsive. On a uniform 36 s mesh, each burn is smeared over whole intervals at low thrust (peak 1.66 N), so the burn shapes are mesh artifacts.
  - The earlier mesh-refinement version (commit `e963d5e`) resolved them, with saturated 10 N burns lasting a few seconds.
  - For the document's initial state and T_max = 10 N, it also reached a lower impulse: J = 131.57 N s at N = 770 nonuniform intervals.
- **Local optimum only.** With the document's initial state and T_max = 10 N, cold starts on different N reach different optima: J = 132.78 at N = 200, 136.80 at N = 300 and 138.08 at N = 400.
- **IPOPT setting.** `bound_relax_factor = 0` is required. With IPOPT's default, `honor_original_bounds` clips the returned thrusts by ~1e-8·T_max. That leaves an ω bias of about 1e-9 s⁻¹ per interval, which grows to about 1e-4 rad of terminal attitude error (measured).

## Outputs (`results/`)

- **Figures** (written by `main` and left open in MATLAB; saved as PNG)
  - `tracking_errors.png`: θ, ω, r, ṙ, with tolerance bands and magnified arrival views.
  - `control_wrench.png`: body-frame force and torque.
  - Thrust: `thrust_histories.png` (per pod), `thrust_burn_zoom.png`, `thrust_heatmap.png`.
  - `impulse.png`: cumulative impulse and per-thruster impulse.
  - `propagation_discrepancy.png`.
  - Static trajectories: `trajectory_formation.png`, `trajectory_SER.png`, and the halo orbit alone in `halo_orbit.png`.
- **Videos** (written by `makeVideo`, drawn live in a MATLAB figure; off-screen in `-batch` mode)
  - `video_SER.mp4` and `video_formation.mp4` share the same frame times.
  - Both show the target formation pyramid (80 × 480 m, apex at 960 m). The view is zoomed on the spacecraft trajectory and the target point, so the pyramid edges are clipped at the box and show where it continues.
  - 1280×720 at 15 fps, 30 s each, at a constant playback speed (×480). The plumes show the peak thrust since the previous frame, so burns shorter than one frame step (32 s) are not skipped.
  - Next to the 3D view, three time histories with a moving cursor, as in the Malladi animation:
    - the eigenaxis attitude error θ_i;
    - the relative position r^f_pi;
    - the thrust magnitude of each of the 20 thrusters (one colour each), drawn as one vertical line per control interval and overlaid, with T_max. Thrusters with the largest total impulse are drawn first, so weaker ones remain visible on top.
  - Every frame shows the elapsed time and the playback rate.
  - The spacecraft follows Fig. 3: upper and lower cylinders, and a central ring carrying the four pods on ±x^i and ±y^i. All 20 nozzles are drawn at r^i_itk and point along −d_k.
  - Each firing thruster gets its own narrow plume. The plume's brightness is ∝ T/T_max, and its length runs from 30% (so any firing thruster stays visible) up to full length at T_max.
  - In the main view the spacecraft is drawn ×3. A true-scale inset in the body frame (top left) shows which nozzles fire.
- **Data**
  - `solution.mat` (v7.3) holds the solution, the propagation, the reference and the parameters.
  - `summary.json`, `settings.json`.
  - CSV files: thrust history, propagated and collocation states.

## Layout

- `src/lifeParameters.m`: every constant of Sec. 4.3, the thruster geometry, solver settings and gates.
- `src/model/`: frames, dynamics (relative, inertial-path check, absolute), state mappings.
- `src/ref/`: halo propagation and Chebyshev fit.
- `src/ocp/`: collocation NLP, initial guess, solve.
- `src/sim/`, `src/post/`: propagation, burns, validation, saving.
- `src/plot/`: figures and videos.
- `tests/FormationTest.m`: formation construction and the initial condition (see above).
- `tests/ModelTest.m`: constants, thruster positive span, frame conventions, halo period, gravity difference against vpa, eqs. (19)–(22) against the inertial path, CasADi against double, terminal equilibrium, quaternion against DCM kinematics, absolute against relative propagation.
