# LIFEguidance

MATLAB implementation of the three single-spacecraft guidance problems at L2
of the LIFE guidance document (Secs. 4–5). One LIFE collector with 20
fixed, throttleable thrusters must reach a point fixed on a rotating virtual
formation near Sun–Earth L2 within t_f = 4 h, with minimum total impulse and a
free arrival time τ.

- **Problem 1** (weakly constrained, fixed phase): only the thrust and arrival-time bounds.
- **Problem 2** (attitude constrained, fixed phase): Problem 1 plus the Sun-angle path constraint (48). The body axis z^i must stay 50°–85° from the anti-Sun direction (θ_i ∈ [95°, 130°] from the Sun).
- **Problem 3** (attitude constrained, free phase): Problem 2 plus the arrival phase φτ of the formation's rotation about the line of sight as a decision variable, eqs. (49)–(56).

All three problems are stated in the non-spinning pointing frame O\* of the
target formation (Sec. 4.1), so the formation phase enters only the terminal
condition (39). Equation numbers refer to the current document (`LIFE_Documents`,
November 4, 2026). The code structure follows the earlier Malladi et al. (2019)
reproduction on branch `malladi-et-al-claude`.

## Usage

Requires MATLAB (R2026a tested; Optimization and Symbolic toolboxes are needed
for the tests only) and an internet connection on first run. CasADi 3.8.1 with
IPOPT is downloaded into `external/` (git-ignored) by `setupCasadi.m`.

```matlab
main                 % solve Problems 1-3, validate, save, compare (~5 min)
phaseSweep           % Problem 3 phase sweep and reduction check, after main (~5 min)
makeVideo            % render the comparison video from the saved runs (~1.5 min)
runtests('tests')    % model consistency checks (~5 s)
```

`problems = [1 2 3]` at the top of `main.m` selects what is solved, and `p.ocp.p3Starts` sets the Problem 3 multi-start phases. Each problem
writes `results/problemK/`. The comparison reads every saved `problemK`, so a
single problem can be re-run on its own. Every problem starts from the same
initial guess (Sec. 5.2) and there is no warm start, so the solutions compare
1:1.

## Method

- **State and model** (`src/model/`)
  - State in the pointing frame O\*, eq. (34): x = [q\*i; ω\*\*i; r\*\*i; ṙ\*\*i], a Hamilton scalar-first quaternion for R\*_i and the position from the formation centre. O\* points at the target star like the formation frame but does not spin (R^c_\* constant, ω\*_c\* = 0).
  - Dynamics (36)–(38): `relativeDynamics.m` implements the formation-frame model (10), (13), (19)–(22) on the CR3BP of eqs. (14)–(15), with cannonball SRP on the spacecraft only and the thruster wrench (29)–(31). Called with `pointingParams(p)` (no spin, r_fp = 0) it is exactly the pointing-frame model; `tests/ModelTest.m` checks both against a separately written absolute model.
  - The dynamics, the initial state x_i,0 and the Sun-angle constraint are independent of the formation phase. The phase enters only the target state x^tgt(φτ), eq. (39) (`targetState.m`).
  - Target-relative quantities (tracking errors, formation-frame plots and videos) come from the solution via R\*_f = R_z(φ(t)) (`toFormation.m`).
  - The gravity difference uses Battin's f(q) form, so no cancellation occurs.
  - The model is type-generic: the same code is used for the NLP (CasADi) and for the propagation (double).
- **Transcription** (`src/ocp/`)
  - One uniform mesh of `p.ocp.N` = 400 intervals.
  - Radau collocation with 3 points per interval and piecewise-constant thrust, solved with IPOPT/MUMPS.
  - The free arrival time enters via s = t/τ, and the objective is eq. (69).
  - Terminal condition (39): ω, r and ṙ equal x^tgt(φτ), and the attitude is imposed as q_x = q_y = 0, q_z cos(φτ/2) − q_w sin(φτ/2) = 0 (Sec. 5.2, invariant to |q|, admits ±q). Problems 1 and 2 substitute φτ = φ̄0 + φ̇τ.
  - Initial guess (Sec. 5.2): τ = t_f, cubic Hermite position and velocity, slerp attitude and linear rate from x_i,0 to x^tgt at the guessed arrival phase, and T = 0.01 T_max. Scaling and tolerances are as in Sec. 5.2.
  - In the NLP, the reference r_bf(t) is a degree-14 Chebyshev fit of the CR3BP halo on [0, t_f] (fit error 6e-5 m). It enters only through the tidal and SRP terms.
- **Validation** (`src/sim/`, `src/post/`)
  - Independent ode89 propagation of the optimized controls.
    - The pointing-frame state is co-integrated with the reference ODE (14), not the fit.
  - Arrival checks (Sec. 5.2): deviations of r, ṙ and ω from x^tgt(φτ) in the pointing frame, and the eigenaxis angle of R^f_i = R_z(φτ)ᵀ R\*_i.
    - RelTol 1e-11 and AbsTol 1e-12 apply to the scaled states, and the integrator restarts at every control switch.
  - The same controls are also propagated through a separately written absolute model (spacecraft and formation origin in SER, inertial attitude). Its gate is 1% of the arrival tolerances.
  - All reported errors, figures and videos are based on the independent propagation.
- **Conventions**
  - R_x, R_y, R_z are active rotations. This is the only reading under which eq. (65) gives ω^b_cb = +Ω e₃ and eq. (66) gives z\* = [cos 65°, sin 65°, 0].
  - Eq. (38) needs the SER-frame velocity difference, which is reconstructed as R^b_\*(ṙ\*\*i + ω\*_b\* × r\*\*i).

## Problem 2: Sun-angle constraint

- **Constraint (47)–(48)**: σ_i = r̂_isᵀ R^b_\* R\*_i e₃ ∈ [cos 130°, cos 95°], see `sunCosine.m` (type-generic). It depends only on the spacecraft's own attitude and the Sun direction, not on the formation phase.
- **NLP**: `buildOcp(p, ref, 2)` appends one inequality row per Radau collocation point. The last point of each interval is the next node, so every node but t = 0 is covered. The Sun direction is taken from the formation centre (Sec. 5.2).
- **Before solving**: `checkSunFeasibility` checks the initial attitude and the formation attitude over [0, t_f].
- **Verification** (`sunAngle.m`, run for both problems, so Problem 1's violation is measured the same way):
  - θ_i along the independent propagation, using the exact Sun direction from the spacecraft;
  - the maximum violation, accepted for Problem 2 if ≤ 0.01°;
  - the time and intervals outside the band by more than that gate;
  - the violation at the collocation points;
  - the error of the formation-centre direction (≈ 7e-10 rad).
- **Active constraint**: a collocation point within 0.01° of a bound. The report lists each one with its time and IPOPT multiplier. Depending on the local optimum, the solution touches a bound at isolated points or rides it along short arcs.

## Problem 3: free arrival phase

- **Phase.** φτ/π is the last decision variable, with bounds [−1, 1] (eq. (56)). It enters only the 12 terminal rows (39), a single sparse column of the constraint Jacobian. Dynamics, initial condition and Sun-angle constraint are those of Problem 2.
- **Post-processing.** After solving, `main` passes `p` with the solution's formation phase at t = 0, φτ − φ̇τ, to validation, saving and plotting. The propagation runs in the phase-independent pointing frame; the tracking errors, the formation-frame figures and the videos are then expressed in the solution's own formation frame.
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

`initialRelativeState.m` returns this state in the pointing frame O\* of the target formation; it does not depend on the target phase. It replaces the fixed initial state (67)–(68) of the document (see Deviations).

**Scenario.** Both sky directions lie 80° from the anti-Sun direction (+x_c at t = 0), on either side of it, so the slew is 160°. Both ends are inside the Problem 2 field of regard, and Problem 1's unconstrained optimum swings z^i out to 146° from anti-Sun. The earlier ±e₂ setup sits at 90°, outside the band at both ends, so Problem 2 would be infeasible there. A ±65° slew was also tried: Problem 1 stays inside the band without the constraint, so it would not show the difference. For the defaults, the initial attitude error is 160°, and the start lies 173 m from the target corner. To change the setup, edit the sky direction, baselines, apex height, spin rate or phase, or the point, then call `deriveFormations(p)`.

`tests/FormationTest.m` checks the construction. The central test maps the initial state into the absolute model and compares attitude, body rate, position and inertial velocity with the initial pyramid. It does this for all 5 points, for the default and for a randomly oriented, non-default initial formation. The ±e₂ setup is kept in the tests (`e2Setup`) to pin the earlier document's eq. (53) and numbers; `skyBasisIsEq66` pins the current eq. (66).

## Results (N = 400, T_max = 3 N, `results/comparison/comparison.json`)

T_max is set to 3 N in `lifeParameters.m`. All problems start from the same initial guess (Sec. 5.2). Problem 3 is the best of its multi-start.

| Quantity | Problem 1 | Problem 2 | Problem 3 |
|---|---|---|---|
| Solver | `Solve_Succeeded`, 157 it., 13 s | `Solve_Succeeded`, 149 it., 16 s | `Solve_Succeeded`, 118 it., 12 s (best start); 123 s for all 8 starts |
| Arrival phase φτ | 60° (prescribed) | 60° (prescribed) | **−5.46°** (free) |
| Arrival time τ | 14400 s, t_f bound active (multiplier 0.50) | 14400 s, t_f bound active (multiplier 0.64) | 14400 s, t_f bound active (multiplier 0.85) |
| Total impulse J | 166.368 N s (Δm = 0.0771 kg) | 174.858 N s (Δm = 0.0810 kg), +5.1% vs P1 | **154.321 N s** (Δm = 0.0715 kg), −11.7% vs P2, −7.2% vs P1 |
| Max thrust | 1.67 N (bound not active) | 1.61 N (bound not active) | 1.88 N (bound not active) |
| Anti-Sun angle of z^i | 30.3°–95.4° | 50.0°–85.0° | 59.5°–85.0° |
| Max Sun violation (propagated) | 19.7° (not imposed) | 1.3e-4° (gate 0.01°, pass) | 3.4e-4° (pass) |
| Time outside the band (> 0.01°) | 8194 s: below 50° at 463–6381 s (down to 30.3°), above 85° at 11664–13932 s (up to 95.4°) | 0 s | 0 s |
| Sun constraint active | – | 23 collocation points on two short arcs: 50° at 6156–6342 s, 85° at 13235–13290 s | 6 points, 85° at 13572–13631 s |
| Way around the Sun line (`sky_plot.png`) | south, cutting inside the band towards anti-Sun | south, along the 50° edge | north, inside the band |
| Terminal ‖Δr‖ / ‖Δṙ‖ | 4.1e-4 m / 5.8e-7 m/s | 1.5e-6 m / 2.3e-10 m/s | 5.9e-7 m / 8.5e-10 m/s |
| Terminal angle ξ / ‖Δω‖ | 3.6e-5 rad / 1.7e-9 s⁻¹ | 1.1e-7 rad (atan2 form) / 1.3e-11 s⁻¹ | 4.7e-8 rad (atan2 form) / 1.3e-11 s⁻¹ |
| Absolute-model cross-check | ≤ 2.0e-8 m, 8.7e-14 rad (pass) | ≤ 1.6e-8 m, 9.5e-14 rad (pass) | ≤ 8.4e-9 m, 3.2e-14 rad (pass) |
| Accepted | yes | yes (including the Sun-angle gate) | yes (including the Sun-angle gate) |

All gates are 1e-2 m, 1e-5 m/s, 1e-4 rad and 1e-7 s⁻¹. All checks of Sec. 4.5 hold: J1 ≤ J2, J3 ≤ J2, Problem 1 violates the Sun-angle constraint, and Problems 2 and 3 satisfy it. With the free phase, the target point arrives elsewhere on its circle about z_f (`trajectory_comparison.png`), which saves more than the constraint costs: J3 < J1.

### Check of the pointing-frame formulation

The pointing-frame formulation replaced a formation-frame one (commit `e135924`), in which the phase entered the dynamics, the Sun constraint and the initial state. Both describe the same physics, so a thrust history must give the same motion in both. This was checked in both directions:

| | Problem 1 | Problem 2 | Problem 3 |
|---|---|---|---|
| Earlier formation-frame optimum, J | 168.401 N s | 174.896 N s | 154.321 N s at −5.46° |
| … validated by the pointing-frame code | accepted, Sun violation 61.0° (as before) | accepted | accepted |
| Pointing-frame optimum (this table), J | 166.368 N s | 174.858 N s | 154.321 N s at −5.46° |
| … validated by the formation-frame code | accepted, Sun violation 19.7° | accepted | accepted |

Problem 3 reaches the same optimum. Problems 1 and 2 reach different, slightly better local optima, because the common initial guess is now interpolated in the pointing frame instead of the spinning formation frame (Sec. 5.2). `tests/ModelTest.m` also propagates one random thrust history in the pointing-frame, formation-frame and absolute models and requires agreement to 1% of the arrival gates.

### Problem 3 multi-start and phase sweep

| Start φτ | Converged φτ | J (N s) |
|---|---|---|
| −180° | −150.1° | 224.16 |
| −135° | **−5.46°** | **154.32** |
| −90°, −45° | −5.15° | 154.42 |
| 0°, 45°, 90°, 135° | 8.2°–8.3° | 156.00–156.03 |

All 8 starts converge (in the formation-frame formulation the −180° start diverged).

**Phase sweep** (`phaseSweep.m`, 4.8 min, `phase_sweep.png`):
- With φτ fixed, the cost has one smooth basin. The minimum is 155.04 N s at 0°; at −22.5° it is 160.85 and at 22.5° it is 158.94.
- J3 = 154.32 N s at −5.46° lies inside that basin and below the sweep minimum, so the multi-start optimum is consistent.
- The worst arrival phases (around ±180°) cost 220–235 N s.

**Reduction check:**
- With φτ fixed at the prescribed 60°, Problem 3 gives J = 174.918 N s with τ = t_f, so it is feasible for Problem 2. Problem 2 gives 174.858 N s.
- The 0.06 N s (0.03%) difference is between local optima of Problem 2: the two NLPs coincide at τ = t_f but are different problems for IPOPT (Problem 2's terminal condition depends on τ), so their iterates end in different nearby optima. Earlier runs found 174.863 and 174.896 N s for the same problem.

Choice of N (measured with the document's initial state and T_max = 10 N, earlier formulation):
with N = 200, all checks pass except the terminal attitude check, which gives
ξ = 2.3e-4 rad. No mesh-convergence study is made (see Deviations).

### Known limitations

- **Burns are not resolved.** The optimum is near-impulsive. On a uniform 36 s mesh, each burn is smeared over whole intervals at low thrust (peak 1.6–1.9 N), so the burn shapes are mesh artifacts.
  - The earlier mesh-refinement version (commit `e963d5e`) resolved them, with saturated 10 N burns lasting a few seconds.
  - For the document's initial state and T_max = 10 N, it also reached a lower impulse: J = 131.57 N s at N = 770 nonuniform intervals.
- **Problems 1 and 2 have several local optima.** They lie within about 1% (Problem 1: 166.37 and 168.40 N s) and 0.2% (Problem 2: 174.60–174.92 N s) of each other, some on opposite sides of the Sun line (the disconnected path sets of Sec. 4.5). Which one IPOPT reaches depends on the initial guess, as the formulation check and the reduction check above show. The 1:1 comparison from the common initial guess therefore does not necessarily give the best optimum of each problem.
- **Problem 3 multi-start.** Each start is capped at `p.ocp.p3MaxWallTime` = 120 s (in the formation-frame formulation one start diverged for more than 40 min); no start reaches the cap now. The multi-start makes `main` take about 5 min; reduce `p.ocp.p3Starts` for faster iteration.
- **Local optimum only.** With the document's initial state and T_max = 10 N, cold starts on different N reach different optima: J = 132.78 at N = 200, 136.80 at N = 300 and 138.08 at N = 400 (earlier formulation).
- **IPOPT setting.** `bound_relax_factor = 0` is required. With IPOPT's default, `honor_original_bounds` clips the returned thrusts by ~1e-8·T_max. That leaves an ω bias of about 1e-9 s⁻¹ per interval, which grows to about 1e-4 rad of terminal attitude error (measured).

### Deviations from the document

- **Scenario (§5.1).** The ±80° slew, the formation-based initial state and the target point at a pyramid corner (above) are kept instead of the document's target pointing (66), initial state (67)–(68) and r_fp = [50, 0, 0] m. The formulation does not depend on this choice; `initialRelativeState` returns the initial state in the pointing frame either way.
- **Mesh (§5.2).** One uniform mesh with N = 400 instead of N = 200 refined until doubling N changes J by less than 1%; no mesh-convergence output.
- **Reference.** A degree-14 Chebyshev fit of the halo instead of a B-spline interpolant (fit error 6e-5 m).
- **Scaling.** Thrust is scaled by T_max = 3 N instead of 10 N; the other scales are those of §5.2.

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
  - `tracking_errors.png`: ξ, ω^f_fi, r^f_pi, ṙ^f_pi (target-relative, formation frame), with tolerance bands and magnified arrival views.
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
  - CSV files: thrust history, propagated states and collocation states. Columns ending `_pt` are the pointing-frame state (34); the propagated file also has the target-relative formation-frame state (`_fp`: q_fi, ω^f_fi, r^f_pi, ṙ^f_pi) and the Sun angle.

## Layout

- `src/lifeParameters.m`: every constant of Sec. 5.1, the thruster geometry, the Sun-angle bounds, solver settings and gates.
- `src/model/`: frames, dynamics (relative, inertial-path check, absolute), the pointing-frame parameters and target state, state mappings.
- `src/ref/`: halo propagation and Chebyshev fit.
- `src/ocp/`: collocation NLP (Problems 1–3), initial guess, solve, Problem 3 multi-start.
- `src/sim/`, `src/post/`: propagation, burns, validation, the Sun angle, saving, the comparison of problems.
- `src/plot/`: figures and videos.
- `tests/FormationTest.m`: formation construction and the initial condition (see above).
- `tests/ModelTest.m`: constants, thruster positive span, frame conventions, halo period, gravity difference against vpa, eqs. (19)–(22) against the inertial path, CasADi against double, the Sun cosine (44), terminal equilibrium, quaternion against DCM kinematics, absolute against relative propagation.
