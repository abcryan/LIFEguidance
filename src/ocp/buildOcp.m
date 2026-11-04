function ocp = buildOcp(p, ref, problem)
%BUILDOCP Minimum-impulse transfer by Radau collocation on p.ocp.N uniform
%   intervals of the normalized time s = t/tau. problem = 1 (default):
%   Problem 1, eqs. (40)-(45); problem = 2: Problem 2, which adds the
%   Sun-angle constraint (48) at every collocation point; problem = 3:
%   Problem 3, eqs. (49)-(56), which adds the arrival phase phiTau.
%   State in the non-spinning pointing frame O* (Sec. 4.1): the dynamics,
%   the initial state and the Sun-angle constraint are phase-independent;
%   phiTau enters only the terminal condition (39).
%   Decision vector w = [tau/tf; X(:); XC(:); U(:)] (Problem 3: and
%   phiTau/pi last) with scaled states X (13 x N+1 nodes), XC (13 x N*d
%   collocation points) and piecewise-constant T/Tmax U (20 x N). Each
%   interval lasts h = tau/N and the objective is (69), h*sum_n sum_k T_k,n.
import casadi.*
if nargin < 3, problem = 1; end
N = p.ocp.N; d = p.ocp.d; nx = 13; nu = p.nThr;
sigma = ones(1, N)/N;
[cRoot, C, D] = radauCoefficients(d);
S = p.ocp.xScale;

% Scaled dynamics dxhat/dt = S^-1 f(t, S xhat, Tmax uhat), eqs. (36)-(38)
pP = pointingParams(p);
t  = SX.sym('t');
xh = SX.sym('x', nx);
uh = SX.sym('u', nu);
f  = relativeDynamics(t, S.*xh, p.Tmax*uh, chebEvalRef(ref, t), pP);
F  = Function('F', {t, xh, uh}, {f./S});
Fm = F.map(N*d);

% Decision variables
nW = 1 + nx*(N + 1) + nx*N*d + nu*N + (problem >= 3);
w = MX.sym('w', nW);
ix.tau = 1;
ix.X   = 1 + (1:nx*(N + 1));
ix.XC  = ix.X(end) + (1:nx*N*d);
ix.U   = ix.XC(end) + (1:nu*N);
ix.phi = [];
if problem >= 3, ix.phi = nW; end
tauH = w(ix.tau);
X  = reshape(w(ix.X), nx, N + 1);
XC = reshape(w(ix.XC), nx, N*d);
U  = reshape(w(ix.U), nu, N);

tau   = p.tf*tauH;
if problem >= 3
    phiTau = pi*w(ix.phi);                    % decision variable, (56)
else
    phiTau = p.target.phase0 + p.target.rate*tau;   % prescribed, (43)
end
sNode = [0, cumsum(sigma)];
sColl = reshape(sNode(1:N) + cRoot(2:end)'*sigma, 1, []);  % column (n-1)*d+j
colInt = kron(1:N, ones(1, d));          % interval of each collocation column
h = tau*sigma;                           % 1xN, MX
Fval = Fm(tau*sColl, XC, U(:, colInt));

% Collocation and continuity
defect = MX(nx, N*d);
xEnd = D(1)*X(:, 1:N);
for j = 1:d
    xp = C(1, j+1)*X(:, 1:N);
    for r = 1:d
        xp = xp + C(r+1, j+1)*XC(:, r:d:end);
    end
    cols = j:d:N*d;
    defect(:, cols) = xp - repmat(h, nx, 1).*Fval(:, cols);
    xEnd = xEnd + D(j+1)*XC(:, j:d:end);
end
cont = X(:, 2:N+1) - xEnd;

% Boundary conditions (42)-(43): x(0) = x_i,0 and x(tau) = x^tgt(phiTau),
% eq. (39); the attitude as q_x = q_y = 0, q_z cos(phiTau/2) -
% q_w sin(phiTau/2) = 0 (Sec. 5.2, invariant to |q|, admits +-q)
x0h = initialRelativeState(p)./S;
xT = S.*X(:, N+1);
xTgt = targetState(phiTau, p);
gT = [xT(2:3); xT(4)*cos(phiTau/2) - xT(1)*sin(phiTau/2); (xT(5:13) - xTgt(5:13))./S(5:13)];
g = [defect(:); cont(:); X(:,1) - x0h; gT];
lbg = zeros(size(g)); ubg = zeros(size(g));

% Sun-angle constraint (48) at the collocation points (the last one of each
% interval is the next node; the initial node is checked before solving),
% with the Sun direction from the formation centre (Sec. 5.2)
ocp.iSun = [];
if problem >= 2
    Sig = Function('Sig', {t, xh}, {sunCosine(t, xh(1:4), chebEvalRef(ref, t), pP)});
    Sm = Sig.map(N*d);
    sig = Sm(tau*sColl, XC);
    ocp.iSun = numel(g) + (1:N*d)';
    g = [g; sig(:)];
    lbg = [lbg; cos(p.sun.thetaMax)*ones(N*d, 1)];
    ubg = [ubg; cos(p.sun.thetaMin)*ones(N*d, 1)];
end

J = p.Tmax*sum(h.*sum(U, 1))/p.ocp.Jscale;

% Bounds: (45) with the numerical lower bound, (44)
lbx = -inf(nW, 1); ubx = inf(nW, 1);
lbx(ix.tau) = p.tauMin/p.tf; ubx(ix.tau) = 1;
lbx(ix.U) = 0; ubx(ix.U) = 1;
lbx(ix.phi) = -1; ubx(ix.phi) = 1;          % (56), phiTau/pi

opts = struct('expand', true, 'print_time', false, 'ipopt', p.ocp.ipopt);
if problem >= 3, opts.ipopt.max_wall_time = p.ocp.p3MaxWallTime; end
ocp.solver = nlpsol('solver', 'ipopt', struct('x', w, 'f', J, 'g', g), opts);
ocp.lbx = lbx; ocp.ubx = ubx;
ocp.lbg = lbg; ocp.ubg = ubg; ocp.problem = problem;
ocp.ix = ix; ocp.N = N; ocp.d = d; ocp.nx = nx; ocp.nu = nu;
ocp.sigma = sigma; ocp.sNode = sNode; ocp.sColl = sColl; ocp.cRoot = cRoot;
ocp.C = C; ocp.D = D; ocp.S = S; ocp.F = F;
ocp.nW = nW; ocp.nG = numel(g);
end
