function ctrl = buildNmpc(p)
%BUILDNMPC NMPC policy (22) as a CasADi/IPOPT problem.
%   Piecewise-constant inputs over N intervals, multiple shooting with RK4
%   (the paper uses two-point collocation; both are direct transcriptions).
%   ctrl.solve(z0) returns the first input and solver info.
import casadi.*

N = p.mpc.N; Ts = p.mpc.Ts; M = p.mpc.M; w = p.w;
nz = 39; nu = 24; nl = size(p.los.A, 1);

% Continuous model and stage cost (20)
z = SX.sym('z', nz);
u = SX.sym('u', nu);
Rct  = reshape(z(1:9), 3, 3);
wct  = z(10:12);
rr   = z(13:18);
xi   = z(19:21);
tauc = p.Btau*u;
ell = w.QR*trace(w.D - w.D*Rct) + wct'*w.Qw*wct + rr'*w.Qrho*rr ...
      + xi'*w.Qxi*xi + tauc'*w.Wtau*tauc + u'*w.Wf*u;
f = Function('f', {z, u}, {relativeDynamics(z, u, p), ell});

% One interval: RK4 on the state and on the cost integral
h = Ts/M;
zk = z; Jk = 0;
for m = 1:M
    [k1, l1] = f(zk, u);
    [k2, l2] = f(zk + h/2*k1, u);
    [k3, l3] = f(zk + h/2*k2, u);
    [k4, l4] = f(zk + h*k3, u);
    zk = zk + h/6*(k1 + 2*k2 + 2*k3 + k4);
    Jk = Jk + h/6*(l1 + 2*l2 + 2*l3 + l4);
end
F = Function('F', {z, u}, {zk, Jk});
ctrl.F = F;

% Optimal control problem
opti = Opti();
Z  = opti.variable(nz, N+1);
U  = opti.variable(nu, N);
S  = opti.variable(nl, N);              % LOS slacks
Z0 = opti.parameter(nz);
slackScale = opti.parameter();          % homotopy factor for cold starts

J = 0;
opti.subject_to(Z(:,1) == Z0);
for k = 1:N
    [zNext, Jint] = F(Z(:,k), U(:,k));
    opti.subject_to(Z(:,k+1) == zNext);
    J = J + Jint;

    % Safeguard keeping iterates where the RK4 grid is accurate; inactive
    % at the solutions (peak |w_c/t| is ~0.02 rad/s)
    opti.subject_to(Z(10:12,k+1) <= p.mpc.wMax);
    opti.subject_to(Z(10:12,k+1) >= -p.mpc.wMax);

    % Soft line-of-sight constraint (18)
    opti.subject_to(p.los.A*Z(13:15,k+1) - p.los.b <= S(:,k));
    opti.subject_to(S(:,k) >= 0);
    J = J + slackScale*Ts*(w.slackLin*sum(S(:,k)) + w.slackQuad*sumsqr(S(:,k)));

    % Thruster gimbal polyhedron (16) and magnitude (17)
    Uk = reshape(U(:,k), 3, 8);
    opti.subject_to(vec(p.thr.A*Uk) <= 0);
    opti.subject_to(sum(Uk.^2, 1)' <= p.thr.fmax^2);
end
rrN = Z(13:18, N+1);
J = J + rrN'*w.P*rrN;                   % terminal cost (21)
opti.minimize(J);

ipoptOpts = struct('print_level', 0, 'max_iter', 1000, 'tol', 1e-6, ...
                   'acceptable_tol', 1e-4);
opti.solver('ipopt', struct('expand', true, 'print_time', false), ipoptOpts);

ctrl.opti = opti;
ctrl.Z = Z; ctrl.U = U; ctrl.S = S; ctrl.Z0 = Z0; ctrl.slackScale = slackScale;
ctrl.solve = @(z0, guess) solveOcp(ctrl, p, z0, guess);
end

function [u0, info] = solveOcp(ctrl, p, z0, guess)
if isempty(guess)
    % Cold start: small nominal thrust along -j-hat, open-loop rollout,
    % then a homotopy from a weakly to the fully penalized LOS constraint
    N = p.mpc.N;
    guess.U = repmat([0; -0.5; 0], 8, N);
    guess.Z = zeros(39, N+1); guess.Z(:,1) = z0;
    for k = 1:N
        guess.Z(:,k+1) = full(ctrl.F(guess.Z(:,k), guess.U(:,k)));
    end
    guess.S = max(p.los.A*guess.Z(13:15,2:end) - p.los.b, 0) + 1e-3;
    [~, info] = solveOnce(ctrl, z0, guess, 1e-3);
    guess = struct('Z', info.Z, 'U', info.U, 'S', info.S);
    tPre = info.solveTime;
    [u0, info] = solveOnce(ctrl, z0, guess, 1);
    info.solveTime = info.solveTime + tPre;
else
    [u0, info] = solveOnce(ctrl, z0, guess, 1);
end
end

function [u0, info] = solveOnce(ctrl, z0, guess, slackScale)
opti = ctrl.opti;
Z = guess.Z; U = guess.U; S = guess.S;
Z(:,1) = z0;
opti.set_value(ctrl.slackScale, slackScale);
opti.set_value(ctrl.Z0, z0);
opti.set_initial(ctrl.Z, Z);
opti.set_initial(ctrl.U, U);
opti.set_initial(ctrl.S, S);

tic;
try
    sol = opti.solve_limited();
    ok = true;
catch
    sol = opti.debug;
    ok = false;
end
info.solveTime = toc;
stats = opti.stats();
info.status = stats.return_status;
info.success = ok && stats.success;
info.iter = stats.iter_count;
info.cost = full(sol.value(opti.f));
info.Z = full(sol.value(ctrl.Z));
info.U = full(sol.value(ctrl.U));
info.S = full(sol.value(ctrl.S));
if ~ok || any(isnan(info.U(:))) || any(isnan(info.Z(:)))
    % Fall back to the (shifted previous) plan used as initial guess
    info.success = false;
    info.Z = Z; info.U = U; info.S = S;
end
u0 = info.U(:,1);
end
