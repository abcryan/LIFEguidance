function ocp = buildOcp(p, ref, sigma, warm)
%BUILDOCP Minimum-impulse transfer, eqs. (35)-(40), by Radau collocation.
%   sigma: interval fractions of the normalized time s = t/tau (sum 1).
%   warm (optional): use the warm-start IPOPT options p.ocp.ipoptWarm.
%   Decision vector w = [tau/tf; X(:); XC(:); U(:)] with scaled states
%   X (13 x N+1 nodes), XC (13 x N*d collocation points) and piecewise-
%   constant T/Tmax U (20 x N). Interval n lasts h_n = tau*sigma_n and the
%   objective is sum_n h_n sum_k T_k,n, which equals (58) on a uniform mesh.
import casadi.*
sigma = sigma(:)'/sum(sigma);
N = numel(sigma); d = p.ocp.d; nx = 13; nu = p.nThr;
[cRoot, C, D] = radauCoefficients(d);
S = p.ocp.xScale;

% Scaled dynamics dxhat/dt = S^-1 f(t, S xhat, Tmax uhat)
t  = SX.sym('t');
xh = SX.sym('x', nx);
uh = SX.sym('u', nu);
f  = relativeDynamics(t, S.*xh, p.Tmax*uh, chebEvalRef(ref, t), p);
F  = Function('F', {t, xh, uh}, {f./S});
Fm = F.map(N*d);

% Decision variables
nW = 1 + nx*(N + 1) + nx*N*d + nu*N;
w = MX.sym('w', nW);
ix.tau = 1;
ix.X   = 1 + (1:nx*(N + 1));
ix.XC  = ix.X(end) + (1:nx*N*d);
ix.U   = ix.XC(end) + (1:nu*N);
tauH = w(ix.tau);
X  = reshape(w(ix.X), nx, N + 1);
XC = reshape(w(ix.XC), nx, N*d);
U  = reshape(w(ix.U), nu, N);

tau   = p.tf*tauH;
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

% Boundary conditions (37)-(38): q_v = 0, w = r = v = 0 at tau
x0h = initialRelativeState(p)./S;
g = [defect(:); cont(:); X(:,1) - x0h; X(2:nx, N+1)];

J = p.Tmax*sum(h.*sum(U, 1))/p.ocp.Jscale;

% Bounds: (39) with the numerical lower bound, (40)
lbx = -inf(nW, 1); ubx = inf(nW, 1);
lbx(ix.tau) = p.tauMin/p.tf; ubx(ix.tau) = 1;
lbx(ix.U) = 0; ubx(ix.U) = 1;

ipopt = p.ocp.ipopt;
if nargin > 3 && warm                   % warm start: begin near the solution
    for fn = fieldnames(p.ocp.ipoptWarm)'
        ipopt.(fn{1}) = p.ocp.ipoptWarm.(fn{1});
    end
end
opts = struct('expand', true, 'print_time', false, 'ipopt', ipopt);
ocp.solver = nlpsol('solver', 'ipopt', struct('x', w, 'f', J, 'g', g), opts);
ocp.lbx = lbx; ocp.ubx = ubx;
ocp.lbg = zeros(size(g)); ocp.ubg = zeros(size(g));
ocp.ix = ix; ocp.N = N; ocp.d = d; ocp.nx = nx; ocp.nu = nu;
ocp.sigma = sigma; ocp.sNode = sNode; ocp.sColl = sColl; ocp.cRoot = cRoot;
ocp.C = C; ocp.D = D; ocp.S = S; ocp.F = F;
ocp.nW = nW; ocp.nG = numel(g);
end
