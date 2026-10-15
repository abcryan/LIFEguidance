# LIFEguidance

MATLAB reproduction of Malladi, Di Cairano & Weiss, *Nonlinear MPC of Coupled
Rotational-Translational Spacecraft Relative Motion* (ACC 2019). A chaser
spacecraft with 8 gimbaled thrusters rendezvous with a tumbling asteroid's
surface rock under thrust, gimbal and line-of-sight constraints.

## Usage

Requires MATLAB (R2026a tested) and an internet connection on first run.

```matlab
main        % simulate 1800 s closed loop (~2 min), save results/ and Figs. 4-7
makeVideo   % render results/inertialFrame.mp4 from the saved run (~4 min)
```

CasADi/IPOPT is downloaded automatically into `external/` (git-ignored) by
`setupCasadi.m`.

## Layout

- `src/paperParameters.m` – all parameters; assumptions not given in the paper are marked `ASSUMED`
- `src/model/` – prediction (relative) and truth (inertial) dynamics
- `src/control/buildNmpc.m` – NMPC problem (CasADi + IPOPT)
- `src/sim/`, `src/plot/` – closed-loop run, summary, figures, animation

## Deviations from the paper

The paper's `W_tau = 9e13` is treated as a typo (9 is used), integral action
only accumulates within 5 m of the goal (anti-windup), and RK4 multiple
shooting replaces collocation. See `paperParameters.m` for details.
