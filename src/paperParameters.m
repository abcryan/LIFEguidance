function p = paperParameters()
%PAPERPARAMETERS Parameters of Malladi, Di Cairano, Weiss (ACC 2019), Sec. IV.
%   Values not stated in the paper are marked "ASSUMED".

% Central body (Earth) and target orbit: circular, 300 km altitude
p.mu  = 3.986004418e14;                 % [m^3/s^2]
p.Re  = 6378137;                        % [m]
p.alt = 300e3;                          % [m]

% Chaser: 1.5 m x 4 m x 1.5 m cuboid
p.mc = 4000;                            % [kg]
p.Jc = diag([6.0833 1.5 6.0833])*1e3;   % [kg m^2]
p.JcInv = inv(p.Jc);
p.dims = [1.5; 4; 1.5];                 % [m]

% Target asteroid. ASSUMED: isotropic inertia, so the tumble rate is
% constant (consistent with the constant steady-state forces in Fig. 7).
p.Jt = eye(3);
p.JtInv = inv(p.Jt);

% Non-center-of-mass feature points
p.rPt = [1.1404; 3.3462; 5.8907];       % rock on target, in F_t [m]
p.rPc = [0; 0; -1.75];                  % end effector on chaser, in F_c [m]

% Thrusters: one on each corner of the cuboid (8 = corners of two faces).
% F_k has j-hat along r_tk/c; nominal thrust is along -j-hat (through c).
% ASSUMED: i-hat = unit(j-hat x z_c), k-hat = i-hat x j-hat.
corners = [ 1  1  1; -1  1  1; -1 -1  1;  1 -1  1;
            1  1 -1; -1  1 -1; -1 -1 -1;  1 -1 -1]';
p.thr.pos = corners .* (p.dims/2);
p.thr.R = zeros(3,3,8);
for k = 1:8
    jHat = p.thr.pos(:,k)/norm(p.thr.pos(:,k));
    iHat = cross(jHat,[0;0;1]); iHat = iHat/norm(iHat);
    kHat = cross(iHat,jHat);
    p.thr.R(:,:,k) = [iHat jHat kHat];  % R_k/c: F_k components -> F_c
end
p.thr.fmax = 20;                        % [N]
p.thr.halfAngle = pi/4;                 % gimbal polyhedron angle [rad]
tg = tan(p.thr.halfAngle);
% A_k f_tk <= 0: square pyramid about -j-hat, eq. (16) with n = 4
p.thr.A = [ 1 tg  0;
           -1 tg  0;
            0 tg  1;
            0 tg -1];

% Force/torque maps (11), (12): fc = Bf*fps, tauc = Btau*fps
p.Bf   = zeros(3,24);
p.Btau = zeros(3,24);
for k = 1:8
    idx = 3*k-2:3*k;
    p.Bf(:,idx)   = p.thr.R(:,:,k);
    p.Btau(:,idx) = skew(p.thr.pos(:,k))*p.thr.R(:,:,k);
end

% Line-of-sight polyhedron (18): pyramid with apex at the rock, axis +z_t.
% ASSUMED: approach axis +z_t (the end effector points along -z_c).
p.los.halfAngle = pi/4;
tl = tan(p.los.halfAngle);
p.los.A = [ 1  0 -tl;
           -1  0 -tl;
            0  1 -tl;
            0 -1 -tl];
p.los.b = zeros(4,1);

% NMPC: N = 8 intervals, horizon T = 480 s
p.mpc.N  = 8;
p.mpc.T  = 480;
p.mpc.Ts = p.mpc.T/p.mpc.N;             % 60 s
p.mpc.M  = 20;                          % RK4 substeps per interval
p.mpc.wMax = 0.2;                       % safeguard bound on |w_c/t| [rad/s]
% ASSUMED anti-windup: with xi(0) = 0 and the 113 m approach, xi = int(rho)
% reaches ~1e4 m s and Qxi*|xi|^2 then dominates Qrho*|rho|^2, driving the
% feature point to the far side of the target. Integrating only near the
% goal (conditional integration) keeps the offset-free property of (14).
p.mpc.xiRadius = 5;                     % integrate xi only for |rho| <= this [m]
p.mpc.xiMax = inf;                      % optional clamp on |xi_i| [m s]

% Cost weights (20), (21)
w.QR   = 70;
w.D    = eye(3);                        % ASSUMED
w.Qw   = 5*eye(3);
w.Qrho = 0.01*eye(6);                   % on [rho; rhodot]
w.Qxi  = 1e-6*eye(3);
% ASSUMED: the paper's Wtau = 9e13*I3 would prevent any attitude
% correction (cost of 1 N m ~ 1e14 vs. 280 for a 180 deg error), so it is
% read as a typo. This value gives the ~200 s attitude convergence of Fig. 6.
w.Wtau = 9*eye(3);
w.Wf   = eye(24);
w.P    = blkdiag(1e2*eye(3), 1e3*eye(3)); % ASSUMED terminal weight
w.slackLin  = 10;                       % ASSUMED soft-LOS penalties
w.slackQuad = 1;
p.w = w;

% Initial conditions
p.ic.Rte = eye(3);
p.ic.wt  = 0.0046*[1; 1; 1];
p.ic.Rce = diag([-1 -1 1]);
p.ic.wc  = zeros(3,1);
p.ic.rho    = [-1.1394; 113.0847; -7.6204];
p.ic.rhoDot = [ 0.0081;   0.0039;  1.3490];
rt = p.Re + p.alt;
p.ic.rt = [rt; 0; 0];
p.ic.vt = [0; sqrt(p.mu/rt); 0];

% Closed-loop simulation
p.sim.Tend = 1800;
p.sim.dtLog = 1;
end
