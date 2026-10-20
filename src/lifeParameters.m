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

%% Frames, eqs. (52)-(54)
p.obliquity = deg2rad(23.43928);        % epsilon
p.ws   = 2*pi/86400;                    % formation spin [rad/s]
p.phi0 = 0;
p.rfp  = [50; 0; 0];                    % target point in O^f [m]

%% Spacecraft, Sec. 4.3
p.m    = 4000;                          % [kg]
p.J    = 1e4*diag([1.625 1.625 2.075]); % [kg m^2]
p.Jinv = inv(p.J);
p.Tmax = 10;                            % [N]
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

%% Initial conditions and arrival, eqs. (55)-(56)
p.ic.Rfi  = rot2(pi/6);                 % Ry(pi/6)
p.ic.wfi  = 1e-4*[1; -1; 1];            % [1/s]
p.ic.rpi  = [100; -50; 20];             % [m]
p.ic.vpi  = [0.01; -0.005; 0];          % [m/s]
p.tf     = 4*3600;                      % arrival deadline [s]
p.tauMin = 1;                           % numerical lower bound on tau [s]

%% Transcription, Sec. 4.3 "Numerical solution and verification"
p.ocp.N0     = 200;                     % initial uniform intervals
p.ocp.d      = 3;                       % Radau collocation points
p.ocp.uGuess = 0.01;                    % T/Tmax initial guess
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
% Mesh refinement (CHOICE where the document is silent)
p.mesh.burnTol   = 1e-3*p.Tmax;         % sum_k T_k above this = firing [N]
p.mesh.relTolJ   = 0.01;                % "doubling N changes J < 1%"
p.mesh.maxLevels = 8;
p.mesh.maxN      = 2000;
p.mesh.locErrTol = 1e-3;                % per-interval error / gate
% Per-burn convergence between successive levels (reported; features that
% do not converge are labelled mesh-dependent)
p.mesh.burnRelTol  = 0.01;              % inertial force/torque integrals
p.mesh.burnFloorF  = 1e-2;              % [N s]
p.mesh.burnFloorT  = 1e-2;              % [N m s]
p.mesh.burnShapeTol = 0.05;             % duration and peak thrust
% IPOPT options for warm-started levels
p.ocp.ipoptWarm = struct('mu_init', 1e-6, 'bound_push', 1e-8, ...
                         'bound_frac', 1e-8);

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
