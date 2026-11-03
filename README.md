# LIFEguidance

MATLAB implementation of the three single-spacecraft guidance problems at L2
of the LIFE guidance document (Secs. 4–5). One LIFE collector with 20
fixed, throttleable thrusters must reach a point fixed on a rotating virtual
formation near Sun–Earth L2 within t_f = 4 h, with minimum total impulse and a
free arrival time τ.

- **Problem 1** (weakly constrained, fixed phase): only the thrust and arrival-time bounds.
- **Problem 2** (attitude constrained, fixed phase): Problem 1 plus the Sun-angle path constraint (45). The body axis z^i must stay 50°–85° from the anti-Sun direction (θ_i ∈ [95°, 130°] from the Sun).
- **Problem 3** (attitude constrained, free phase): Problem 2 plus the arrival phase φτ of the formation's rotation about the line of sight as a decision variable, eqs. (46)–(55).

The code structure follows the earlier Malladi et al. (2019) reproduction on
branch `malladi-et-al-claude`.

## Usage

Requires MATLAB (R2026a tested; Optimization and Symbolic toolboxes are needed
for the tests only) and an internet connection on first run. CasADi 3.8.1 with
IPOPT is downloaded into `external/` (git-ignored) by `setupCasadi.m`.

```matlab
main                 % solve Problems 1-3, validate, save, compare (~7 min)
phaseSweep           % Problem 3 phase sweep and reduction check, after main (~4.5 min)
makeVideo            % render the comparison video from the saved runs (~1.6 min)
runtests('tests')    % model consistency checks (~5 s)
```

`problems = [1 2 3]` at the top of `main.m` selects what is solved, and `p.ocp.p3Starts` sets the Problem 3 multi-start phases. Each problem
writes `results/problemK/`. The comparison reads every saved `problemK`, so a
single problem can be re-run on its own. Every problem starts from the same
initial guess (Sec. 5.2) and there is no warm start, so the solutions compare
1:1.

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

## Problem 2: Sun-angle constraint

- **Constraint (45)**: σ_i = r̂_isᵀ R^b_f R(q) e₃ ∈ [cos 130°, cos 95°], see `sunCosine.m` (type-generic).
- **NLP**: `buildOcp(p, ref, 2)` appends one inequality row per Radau collocation point. The last point of each interval is the next node, so every node but t = 0 is covered. The Sun direction is taken from the formation centre (Sec. 5.2).
- **Problem 1 is unchanged**: `buildOcp(p, ref, 1)` builds exactly the previous NLP.
- **Before solving**: `checkSunFeasibility` checks the initial attitude and the formation attitude over [0, t_f].
- **Verification** (`sunAngle.m`, run for both problems, so Problem 1's violation is measured the same way):
  - θ_i along the independent propagation, using the exact Sun direction from the spacecraft;
  - the maximum violation, accepted for Problem 2 if ≤ 0.01°;
  - the time and intervals outside the band by more than that gate;
  - the violation at the collocation points;
  - the error of the formation-centre direction (≈ 7e-10 rad).
- **Active constraint**: a collocation point within 0.01° of a bound. The report lists each one with its time and IPOPT multiplier. Depending on the local optimum, the solution touches a bound at isolated points or rides it along short arcs.

## Problem 3: free arrival phase

- **Phase.** Every phase-dependent quantity in the code goes through the target formation's phase at t = 0, `p.target.phase0`.
  - In the NLP it is an argument of the dynamics and Sun-angle functions: fixed for Problems 1 and 2, and φτ − φ̇τ (eq. (46)) for Problem 3.
  - φτ/π is the last decision variable, with bounds [−1, 1] (eq. (55)).
- **Initial condition (47).** The spacecraft's physical initial state is fixed: it is defined absolutely by the initial formation `p.init`. Its relative state for a shifted phase is `shiftPhase`, a closed form that works symbolically. `tests/FormationTest.m` checks it against `initialRelativeState` with the shifted phase.
- **Post-processing.** After solving, `main` passes `p` with the solution's `phase0` to validation, saving and plotting. Propagation, the absolute-model check, the Sun angle and all figures are then evaluated in the solution's own formation frame.
- **Multi-start** (`solveMultiStart`, Sec. 5.2): 8 solves from φτ ∈ {−π + mπ/4}, each from the common initial guess with the boundary states of that phase. The best `Solve_Succeeded` solution is kept, and `results/problem3/multistart.json` lists every start.
- **Phase sweep** (`phaseSweep.m`, run after `main`):
  - 16 solves with φτ fixed, which must not undercut the multi-start optimum;
  - one solve at Problem 2's prescribed arrival phase φ̄0 + φ̇t_f, the reduction check: Problem 3 with this phase is Problem 2;
  - writes `results/problem3/phase_sweep.json` and `phase_sweep.png`.

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
| Sky direction | [cos 80°, sin 80°, 0] | [cos 80°, −sin 80°, 0] |
| Spin | 2π/24 h about its own z | 2π/24 h about its own z |
| Point used | 1 → target point r^f_fp = [40, 240, 0] m | 1 → spacecraft start |

The initial formation shares the target origin at t = 0, in both position and velocity. The spacecraft starts at the chosen point of the initial formation, with its body axes aligned with the initial frame, and co-rotates with it. That means:

- R^c_i(0) = R^c_f0(0);
- ω^i_ci(0) = rate·e₃;
- the inertial velocity relative to the origin is ω_cf0 × r.

This replaces the fixed initial state of eqs. (55)–(56). See `initialRelativeState.m`.

**Scenario.** Both sky directions lie 80° from the anti-Sun direction (+x_c at t = 0), on either side of it, so the slew is 160°. Both ends are inside the Problem 2 field of regard, and Problem 1's unconstrained optimum swings z^i out to 146° from anti-Sun. The earlier ±e₂ setup sits at 90°, outside the band at both ends, so Problem 2 would be infeasible there. A ±65° slew was also tried: Problem 1 stays inside the band without the constraint, so it would not show the difference. For the defaults, the initial attitude error is 160°, and the start lies 173 m from the target corner. To change the setup, edit the sky direction, baselines, apex height, spin rate or phase, or the point, then call `deriveFormations(p)`.

`tests/FormationTest.m` checks the construction. The central test maps the initial state into the absolute model and compares attitude, body rate, position and inertial velocity with the initial pyramid. It does this for all 5 points, for the default and for a randomly oriented, non-default initial formation. The ±e₂ setup is kept in the tests (`e2Setup`) to pin eq. (53) and the earlier numbers.

## Results (N = 400, T_max = 3 N, `results/comparison/comparison.json`)

T_max is set to 3 N in `lifeParameters.m`. All problems start from the same initial guess (Sec. 5.2). Problem 3 is the best of its multi-start.

| Quantity | Problem 1 | Problem 2 | Problem 3 |
|---|---|---|---|
| Solver | `Solve_Succeeded`, 217 it., 20 s | `Solve_Succeeded`, 161 it., 19 s | `Solve_Succeeded`, 119 it., 12 s (best start); 259 s for all 8 starts |
| Arrival phase φτ | 60° (prescribed) | 60° (prescribed) | **−5.46°** (free) |
| Arrival time τ | 14400 s, t_f bound active (multiplier 0.59) | 14400 s, t_f bound active (multiplier 0.65) | 14400 s, t_f bound active (multiplier 0.85) |
| Total impulse J | 168.401 N s (Δm = 0.0781 kg) | 174.896 N s (Δm = 0.0811 kg), +3.9% vs P1 | **154.321 N s** (Δm = 0.0715 kg), −11.8% vs P2, −8.4% vs P1 |
| Max thrust | 1.73 N (bound not active) | 1.66 N (bound not active) | 1.88 N (bound not active) |
| Anti-Sun angle of z^i | 49.7°–146.0° | 50.0°–85.0° | 59.5°–85.0° |
| Max Sun violation (propagated) | 61.0° (not imposed) | < 1e-3° (gate 0.01°, pass) | < 1e-3° (pass) |
| Time outside the band (> 0.01°) | 9783 s, mostly 4104–13842 s | 0 s | 0 s |
| Sun constraint active | – | 22 collocation points on two short arcs: 50° at 6126–6306 s, 85° at 13182–13248 s | 6 points, 85° at 13572–13631 s |
| Way around the Sun line (`sky_plot.png`) | north of anti-Sun, out to 146° | south of anti-Sun | north, inside the band |
| Terminal ‖r‖ / ‖ṙ‖ | 4.4e-4 m / 4.6e-7 m/s | 1.8e-6 m / 4.6e-10 m/s | 6.4e-7 m / 9.2e-10 m/s |
| Terminal angle / ‖ω‖ | 2.4e-5 rad / 3.7e-9 s⁻¹ | 6.4e-8 rad (atan2 form) / 9.9e-12 s⁻¹ | 5.1e-8 rad (atan2 form) / 1.3e-11 s⁻¹ |
| Absolute-model cross-check | ≤ 2.5e-8 m, 1.7e-13 rad (pass) | ≤ 1.1e-8 m, 6.7e-14 rad (pass) | ≤ 9.7e-9 m, 2.2e-13 rad (pass) |
| Accepted | yes | yes (including the Sun-angle gate) | yes (including the Sun-angle gate) |

All gates are 1e-2 m, 1e-5 m/s, 1e-4 rad and 1e-7 s⁻¹. All checks of Sec. 4.5 hold: J1 ≤ J2, J3 ≤ J2, Problem 1 violates the Sun-angle constraint, and Problems 2 and 3 satisfy it. With the free phase, the target point arrives elsewhere on its circle about z_f (`trajectory_comparison.png`), which saves more than the constraint costs: J3 < J1.

### Problem 3 multi-start and phase sweep

| Start φτ | Converged φτ | J (N s) |
|---|---|---|
| −180° | – (stopped at the 120 s cap, diverged in IPOPT's restoration phase) | – |
| −135° | −148.9° | 223.85 |
| −90°, −45° | −5.3°, **−5.5°** | 154.38, **154.32** |
| 0°, 45°, 90°, 135° | 8.2°–8.4° | 156.02–156.03 |

**Phase sweep** (`phaseSweep.m`, 4.5 min, `phase_sweep.png`):
- With φτ fixed, the cost has one smooth basin. The minimum is 155.05 N s at 0°; at −22.5° it is 160.31 and at 22.5° it is 158.70.
- J3 = 154.32 N s at −5.46° lies inside that basin and below the sweep minimum, so the multi-start optimum is consistent.
- The worst arrival phases (±180°) cost 225–230 N s.

**Reduction check:**
- With φτ fixed at the prescribed 60°, Problem 3 gives J = 174.863 N s. Problem 2 gives 174.896 N s.
- That solution has τ = t_f and φ(0) = 3e-9 rad, so it is a Problem 2 solution. Problem 2's NLP warm-started from it converges to the same 174.863 N s, and Problem 2's validation accepts it.
- So the 0.03 N s difference is a second, slightly better local optimum of Problem 2, not an implementation discrepancy.

Choice of N (measured with the document's initial state and T_max = 10 N):
with N = 200, all checks pass except the attitude check (60), which gives
θ = 2.3e-4 rad.

### Known limitations

- **Burns are not resolved.** The optimum is near-impulsive. On a uniform 36 s mesh, each burn is smeared over whole intervals at low thrust (peak 1.66 N), so the burn shapes are mesh artifacts.
  - The earlier mesh-refinement version (commit `e963d5e`) resolved them, with saturated 10 N burns lasting a few seconds.
  - For the document's initial state and T_max = 10 N, it also reached a lower impulse: J = 131.57 N s at N = 770 nonuniform intervals.
- **Problem 2 has several local optima.** Some lie on opposite sides of the Sun line (the disconnected path sets of Sec. 4.5), and one is 0.03 N s below the reported one (see the reduction check above).
  - From the common initial guess (the reported result), z^i passes south of anti-Sun: J = 174.896 N s.
  - Warm-started from the Problem 1 solution, it passes north, as Problem 1 does: J = 174.595 N s. That is 0.3 N s (0.17%) lower, with 6 touch points and a maximum violation of 0.008°.
  - Both are accepted. The 1:1 comparison therefore does not give the best Problem 2 optimum found.
- **Problem 3 multi-start.** One of the 8 starts (−180°) diverges and is stopped by the per-start cap `p.ocp.p3MaxWallTime` = 120 s; uncapped it ran for more than 40 min. Fixing φτ = −180° in the sweep converges normally, so the failure comes from the start, not from that phase. The multi-start makes `main` take about 7 min; reduce `p.ocp.p3Starts` for faster iteration.
- **Local optimum only.** With the document's initial state and T_max = 10 N, cold starts on different N reach different optima: J = 132.78 at N = 200, 136.80 at N = 300 and 138.08 at N = 400.
- **IPOPT setting.** `bound_relax_factor = 0` is required. With IPOPT's default, `honor_original_bounds` clips the returned thrusts by ~1e-8·T_max. That leaves an ω bias of about 1e-9 s⁻¹ per interval, which grows to about 1e-4 rad of terminal attitude error (measured).

## Outputs (`results/`)

- **Comparison** (`results/comparison/`, written by `main`)
  - `comparison.json` and `comparison.csv`, one row per problem: τ*, J, Δm, φτ* (prescribed or free), active bounds, Sun violation, time outside, active points, solve time (Problem 3 also over all starts) and acceptance. The same table is printed together with the consistency checks (Sec. 4.5, and the phase sweep when present).
  - `sun_angle_comparison.png`: the anti-Sun angle against t for every problem. The field of regard is green, Problem 1's violations are shaded red, and the active points are marked.
  - `sky_plot.png`: the direction of z^i about the Sun line. The radius is the anti-Sun angle and ecliptic north is up; the admissible annulus is green.
  - `tracking_comparison.png` and `impulse_comparison.png` (cumulative and per thruster).
  - `trajectory_comparison.png`: the paths in the prescribed formation frame with z^i arrows. Problem 3 is rotated into this frame, so both target points show on their circle about z_f.
  - `phase_sweep.png` (from `phaseSweep`): J against fixed φτ, the multi-start results, φτ* and Problem 2's prescribed phase.
  - `video_comparison.mp4` (`makeVideo`, 30 s, ×480):
    - top left, one attitude view per problem, in SER axes with the same camera for all: across the Sun line, anti-Sun up. Each view shows the spacecraft at true scale with plumes, the Sun direction, the field of regard as a green zone of the sphere, and z^i, green inside and red outside;
    - bottom left, a formation-centred close-up with SER-parallel axes. It shows the prescribed virtual formation (red, Problems 1 and 2) and Problem 3's formation (purple, rotated by its phase), each with its moving target point, and all three paths;
    - right, a live sky plot, info boxes (anti-Sun angle, inside or outside, impulse, arrival phase), the anti-Sun angle and the cumulative impulse with a cursor.
- **Per problem** (`results/problemK/`): the figures and data below. `sun_angle.png` shows θ_i against the band, with the intervals outside it and the active points. `summary.json` holds a `sun` block and the arrival phase. Problem 3 also writes `multistart.json`, and `phaseSweep` adds `phase_sweep.json` and `.png`. The two earlier videos per problem are rendered with `video.perProblem = true` in `makeVideo.m`; each uses its own formation phase.
- **Figures** (written by `main` and left open in MATLAB; saved as PNG)
  - `tracking_errors.png`: θ, ω, r, ṙ, with tolerance bands and magnified arrival views.
  - `control_wrench.png`: body-frame force and torque.
  - Thrust: `thrust_histories.png` (per pod), `thrust_burn_zoom.png`, `thrust_heatmap.png`.
  - `impulse.png`: cumulative impulse and per-thruster impulse.
  - `propagation_discrepancy.png`.
  - Static trajectories: `trajectory_formation.png`, `trajectory_SER.png`, and the halo orbit alone in `halo_orbit.png`.
- **Videos** (written by `makeVideo`, drawn live in a MATLAB figure; off-screen in `-batch` mode)
  - `video_SER.mp4` and `video_formation.mp4` share the same frame times.
  - Both show the target formation pyramid (80 × 480 m, apex at 960 m) in red. For reference only, they also show the initial formation (25 × 150 m, apex 300 m) in grey, fixed at its t = 0 pose for the whole video. The view is zoomed on the spacecraft trajectory and the target point, so the pyramid edges are clipped at the box and show where it continues.
  - 1920×1080 at 30 fps and MPEG-4 quality 95, 30 s each, at a constant playback speed (×480). These settings are at the top of `makeVideo.m`. The plumes show the peak thrust since the previous frame, so burns shorter than one frame step (32 s) are not skipped.
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
  - CSV files: thrust history, propagated states (with the Sun angle) and collocation states.

## Layout

- `src/lifeParameters.m`: every constant of Sec. 5.1, the thruster geometry, the Sun-angle bounds, solver settings and gates.
- `src/model/`: frames, dynamics (relative, inertial-path check, absolute), state mappings.
- `src/ref/`: halo propagation and Chebyshev fit.
- `src/ocp/`: collocation NLP (Problems 1–3), initial guess, solve, Problem 3 multi-start.
- `src/sim/`, `src/post/`: propagation, burns, validation, the Sun angle, saving, the comparison of problems.
- `src/plot/`: figures and videos.
- `tests/FormationTest.m`: formation construction and the initial condition (see above).
- `tests/ModelTest.m`: constants, thruster positive span, frame conventions, halo period, gravity difference against vpa, eqs. (19)–(22) against the inertial path, CasADi against double, the Sun cosine (44), terminal equilibrium, quaternion against DCM kinematics, absolute against relative propagation.
