function p = lifeParameters()
%LIFEPARAMETERS Benchmark of Guidance_Formulation_LIFE.pdf, Sec. 4.3.
%   Equation numbers refer to that document. Choices not fixed by the
%   document are marked "CHOICE".
%
%   Rotation convention: Rx, Ry, Rz (rot1, rot2, rot3 in the code) are
%   active rotations exp(a*[e_k]x).
%   This is the only convention for which eq. (52) gives w^b_cb = +Omega*e3
%   with z_b on the ecliptic pole, and eq. (53) gives R^c_f*e3 = +e2.

%% CR3BP reference, eqs. (47)-(51)
p.dse  = 149597870700;                  % Sun-Earth distance [m]
p.mus  = 1.32712440041279419e20;        % [m^3/s^2]
p.mue  = 3.98600435507e14;              % [m^3/s^2]
p.eta  = p.mue/(p.mus + p.mue);         % (49)
p.Om   = sqrt((p.mus + p.mue)/p.dse^3); % (50) SER rotation rate [1/s]
p.rbs  = [-p.eta*p.dse; 0; 0];          % Sun in SER [m]
p.rbe  = [(1 - p.eta)*p.dse; 0; 0];     % Earth in SER [m]
p.K    = [0 1 0; -1 0 0; 0 0 0];        % coupling matrix of (14)
p.haloIC = [1.008311026993092; 0; 0.000800000000000; ...
            0; 0.009991106037634; 0];   % (51), normalized [r; v]
p.haloPeriodDays = 180.29945;           % stated approximate period

%% Frames, eq. (52), and the virtual formations
p.obliquity = deg2rad(23.43928);        % epsilon
% A virtual formation is a pyramid: four collectors at the corners of a
% nulling (x) x imaging (y) rectangle in its z = 0 plane, centred on the
% frame origin, and the combiner at the apex, height along +z. The frame
% points at the sky direction (unit vector in the inertial frame c) and
% spins about it: R^c = skyBasis(sky)*Rz(phase0 + rate*t). Points 1-4 are
% the corners (+,+), (-,+), (-,-), (+,-) in (x, y), point 5 is the apex.
% Target formation (the spacecraft must reach p.target.point of it). With
% sky = e2, phase0 = 0 and a 24 h spin it is eq. (53).
% CHOICE: the target and initial skies lie 80 deg from the anti-Sun
% direction (+x_c at t = 0) on either side of it, a 160 deg slew. Both ends
% are inside the field of regard of Problem 2 (the default +-e2 sit at
% 90 deg, outside it), and the unconstrained optimum (Problem 1) swings z^i
% to ~146 deg from anti-Sun, so the Sun-angle constraint matters.
p.target = struct('sky', [cosd(80); sind(80); 0], 'phase0', 0, 'rate', 2*pi/86400, ...
                  'nulling', 80, 'imaging', 480, 'height', 960, 'point', 1);
% Initial formation: exists only at t = 0, shares the target origin, and
% defines the initial condition: the spacecraft sits at p.init.point with
% body axes aligned with this frame, co-rotating with it (replaces the
% fixed initial state of eqs. (55)-(56)).
p.init = struct('sky', [cosd(80); -sind(80); 0], 'phase0', 0, 'rate', 2*pi/86400, ...
                'nulling', 25, 'imaging', 150, 'height', 300, 'point', 1);
p = deriveFormations(p);

%% Fleet of Problem 4, Sec. 5.1 "Additions for Problem 4"
% Spacecraft i of the document -> point of pyramidPoints (whose corner
% order is (+,+), (-,+), (-,-), (+,-)): collectors (+,+), (+,-), (-,+),
% (-,-) and the combiner at the apex. The same point in the initial and
% the target formation (fixed assignment). Spacecraft 1 is the one of
% Problems 1-3. override{i}: spacecraft-specific parameters (fields of p,
% e.g. m, J, Asrp, sun), empty: identical spacecraft.
p.fleet.point = [1 4 2 3 5];
p.fleet.name = ["collector (+,+)", "collector (+,-)", "collector (-,+)", ...
                "collector (-,-)", "combiner"];
p.fleet.override = cell(1, 5);

%% Spacecraft, Sec. 4.3
p.m    = 4000;                          % [kg]
p.J    = 1e4*diag([1.625 1.625 2.075]); % [kg m^2]
p.Jinv = inv(p.J);
p.Tmax = 3;                             % [N] (document: 10 N; reduced on request)
p.nThr = 20;
p.thr  = thrusterGeometry();
p.Bf   = p.thr.dir;                     % force map (29)
p.Btau = zeros(3, p.nThr);              % torque map (30)
for k = 1:p.nThr
    p.Btau(:,k) = cross(p.thr.pos(:,k), p.thr.dir(:,k));
end
p.B = [p.Bf; p.Btau];                   % (46)

%% Solar radiation pressure, eq. (57)
p.Ssun = 1361; p.Cr = 1.5; p.Asrp = 20; p.clight = 299792458;
p.beta = p.Ssun*p.Cr*p.Asrp*p.dse^2/(p.clight*p.m);

%% Sun-angle constraint of Problem 2, eqs. (44)-(45)
% The angle theta_i between the body vector n^i and the Sun direction must
% stay in [thetaMin, thetaMax]: the anti-Sun field of regard [50, 85] deg
% (the document's 45 deg is a typo, corrected by the author).
p.sun.n        = [0; 0; 1];             % n^i = e3
p.sun.thetaMin = deg2rad(95);           % anti-Sun angle <= 85 deg
p.sun.thetaMax = deg2rad(130);          % anti-Sun angle >= 50 deg
% CHOICE: the NLP imposes the band this much inside [thetaMin, thetaMax]
% (back-off for the motion between collocation points). Without it, at
% N = 400, d = 5 the combiner (Problem 4) touched the 50 deg bound at two
% collocation points 10 s apart while slewing fast, and the propagated
% angle dipped to 49.982 deg between them, 0.018 deg past the gate.
% The verification and its gate below use the band itself.
p.sun.marginDeg = 0.05;                 % [deg]
% Problem 2 verification: max violation along the propagation (Sec. 5.2)
p.gate.sunDeg  = 0.01;                  % [deg]

%% Arrival
p.tf     = 4*3600;                      % arrival deadline [s]
p.tauMin = 1;                           % numerical lower bound on tau [s]

%% Transcription, Sec. 4.3 "Numerical solution and verification"
% Uniform intervals. CHOICE: the document's initial N = 200 fails the
% attitude check (60) (theta = 2.3e-4 rad). N = 400 with d = 3 Radau
% points passes for one collector (Problems 1-3), but not for the combiner
% of Problem 4: it slews 160 deg while flying 1.25 km, at up to
% 0.012-0.020 rad/s, i.e. 0.4-0.7 rad per 36 s interval, which three
% collocation points do not resolve (measured at the Problem 4 optimum:
% |dr| 1.9e-2 m, xi 6.9e-4 rad, |dw| 1.4e-7 1/s; N = 500 and 600 still
% fail for its faster-turning local optima). CHOICE: d = 5 Radau points
% (order 9 instead of 5) at N = 400: the combiner passes with large
% margins, also its fastest-turning optimum (0.019 rad/s: xi 3.8e-7 rad,
% |dr| 1.5e-5 m), with 36% more variables. One mesh for all problems
% keeps them comparable.
p.ocp.N      = 400;
p.ocp.d      = 5;                       % Radau collocation points
p.ocp.uGuess = 0.01;                    % T/Tmax initial guess
% Problem 3 multi-start arrival phases, Sec. 5.2: -pi + m*pi/4, m = 0..7
p.ocp.p3Starts = -pi + (0:7)*pi/4;
% CHOICE: wall-time cap per Problem 3 solve [s] (a normal start converges
% in 15-30 s; the start at -pi diverged in IPOPT's restoration phase with
% MUMPS reallocating memory, > 40 min uncapped). A capped start counts as
% not converged and is reported; the best converged start is kept.
p.ocp.p3MaxWallTime = 120;
% State scaling [q; w; r; v]: quaternion unscaled, 1e-3 1/s, 100 m, 0.1 m/s
p.ocp.xScale = [1; 1; 1; 1; 1e-3*ones(3,1); 100*ones(3,1); 0.1*ones(3,1)];
p.ocp.Jscale = 100;                     % CHOICE: objective scale [N s]
% bound_relax_factor = 0: IPOPT's default (1e-8) lets T/Tmax reach -1e-8,
% and honor_original_bounds then clips the returned controls to [0, 1].
% The clipped controls no longer match the collocated states: ~1e-7 N per
% thruster is enough to bias w by ~1e-9 1/s per interval, which integrates
% to ~1e-4 rad over the coast (measured). With 0 the iterates stay feasible.
p.ocp.ipopt  = struct('tol', 1e-8, 'constr_viol_tol', 1e-8, ...
                      'max_iter', 3000, 'honor_original_bounds', 'yes', ...
                      'bound_relax_factor', 0, ...
                      'linear_solver', 'mumps', 'print_level', 5, ...
                      'print_frequency_iter', 50);
% Problem 4 (Sec. 5.2): fixed-phase sweep at tau = tf on 16 phases (plus
% Problem 2's prescribed phase, an upper bound), two attitude-class guesses
% per subproblem through an intermediate z^i tilted from anti-Sun towards
% +-z_c. CHOICE: tilted to the middle of the field of regard (67.5 deg from
% anti-Sun); z_c itself lies 90 deg from anti-Sun, outside the band.
p.ocp.p4.sweepPhases = -pi + (0:15)*pi/8;
p.ocp.p4.viaAntiSunDeg = 67.5;
p.ocp.p4.nJointStarts = 3;              % joint NLP from the best grid points,
                                        % primal start (continuationOptions)
% CHOICE: finite-difference check of dJ_i/dphi: slope of a quadratic fit
% through re-solves at phi* + m*fdStepDeg, m = -2..2 [deg] (solveFleet)
p.ocp.p4.fdStepDeg = 0.25;
% CHOICE: wall-time cap per joint solve [s] (measured 91-220 s)
p.ocp.p4MaxWallTime = 900;
% Parallel pool for the Problem 4 sweep and validation (process workers)
p.par.nWorkers = 8;
% A burn is a run of intervals with sum_k T_k above this (CHOICE) [N]
p.burnTol = 1e-3*p.Tmax;

%% Forward propagation and acceptance, eqs. (59)-(60)
p.prop.RelTol = 1e-11;                  % on scaled states
p.prop.AbsTol = 1e-12;
p.gate.r     = 1e-2;                    % [m]
p.gate.v     = 1e-5;                    % [m/s]
p.gate.theta = 1e-4;                    % [rad]
p.gate.w     = 1e-7;                    % [1/s]
% CHOICE: absolute-model cross-check must agree to 1% of the gates above
p.gate.absRel = 0.01;

%% Propellant estimate
p.Isp = 220; p.g0 = 9.80665;
end

function thr = thrusterGeometry()
%THRUSTERGEOMETRY Positions and force directions in O^i, eqs. (41)-(45).
delta = 0.35;
rho = sqrt(3.8^2 - 2*delta^2);          % (43)
st = [1 1; 1 -1; -1 1; -1 -1];          % (s,t) order of the corner thrusters
thr.pos = zeros(3, 20);
thr.dir = zeros(3, 20);
thr.pod = zeros(1, 20);
thr.label = strings(1, 20);
e3 = [0; 0; 1];
k = 0;
for j = 0:3
    alpha = j*pi/2;                     % (41)
    a = [cos(alpha); sin(alpha); 0];    % (42)
    b = [-sin(alpha); cos(alpha); 0];
    k = k + 1;                          % radial thruster (44)
    thr.pos(:,k) = 3.8*a;
    thr.dir(:,k) = -a;
    thr.pod(k) = j + 1;
    thr.label(k) = sprintf('P%d-R', j + 1);
    for c = 1:4                         % corner thrusters (45)
        s = st(c,1); t = st(c,2);
        k = k + 1;
        thr.pos(:,k) = rho*a + s*delta*b + t*delta*e3;
        thr.dir(:,k) = -(a + s*b + t*e3)/sqrt(3);
        thr.pod(k) = j + 1;
        thr.label(k) = sprintf('P%d-%s%s', j + 1, sgn(s), sgn(t));
    end
end
thr.rho = rho;
thr.delta = delta;
end

function c = sgn(x)
if x > 0, c = '+'; else, c = '-'; end
end
