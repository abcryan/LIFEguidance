function blk = ocpBlock(c, p, tau, phiTau, wX, wXC, wU, par)
%OCPBLOCK Constraints and objective of one spacecraft in the collocation
%   NLP (Problems 1-4): Radau collocation and continuity, the initial
%   condition (42), the terminal condition (39) at the arrival phase
%   phiTau and, if c.Sm is set, the Sun-angle constraint (48) at every
%   collocation point. tau and phiTau are MX expressions (shared by all
%   spacecraft in Problem 4); wX, wXC, wU the block's scaled states and
%   T/Tmax; par = [x_i,0./S; r_fp] (16x1) its initial state and target
%   point, numeric or NLP parameters (MX, one compiled solver for every
%   spacecraft with the same physics).
%   Returns g, lbg, ubg, the scaled objective term J (impulse / Jscale) and
%   the local row indices iTerm (the 12 terminal rows) and iSun.
import casadi.*
N = c.N; d = c.d; nx = c.nx; nu = c.nu; S = c.S;
X  = reshape(wX, nx, N + 1);
XC = reshape(wXC, nx, N*d);
U  = reshape(wU, nu, N);
h = tau*c.sigma;                         % 1xN, MX
Fval = c.Fm(tau*c.sColl, XC, U(:, c.colInt));

% Collocation and continuity
defect = MX(nx, N*d);
xEnd = c.D(1)*X(:, 1:N);
for j = 1:d
    xp = c.C(1, j+1)*X(:, 1:N);
    for r = 1:d
        xp = xp + c.C(r+1, j+1)*XC(:, r:d:end);
    end
    cols = j:d:N*d;
    defect(:, cols) = xp - repmat(h, nx, 1).*Fval(:, cols);
    xEnd = xEnd + c.D(j+1)*XC(:, j:d:end);
end
cont = X(:, 2:N+1) - xEnd;

% Boundary conditions: x(0) = x_i,0 and x(tau) = x^tgt(phiTau), eq. (39);
% the attitude as q_x = q_y = 0, q_z cos(phiTau/2) - q_w sin(phiTau/2) = 0
% (Sec. 5.2, invariant to |q|, admits +-q)
x0h = par(1:nx);
pT = p; pT.rfp = par(nx + (1:3));
xT = S.*X(:, N+1);
xTgt = targetState(phiTau, pT);
gT = [xT(2:3); xT(4)*cos(phiTau/2) - xT(1)*sin(phiTau/2); (xT(5:13) - xTgt(5:13))./S(5:13)];
g = [defect(:); cont(:); X(:,1) - x0h];
blk.iTerm = numel(g) + (1:12)';
g = [g; gT];
lbg = zeros(size(g)); ubg = zeros(size(g));

% Sun-angle constraint (48) at the collocation points (the last one of each
% interval is the next node; the initial node is checked before solving),
% with the Sun direction from the formation centre (Sec. 5.2)
blk.iSun = [];
if ~isempty(c.Sm)
    sig = c.Sm(tau*c.sColl, XC);
    blk.iSun = numel(g) + (1:N*d)';
    g = [g; sig(:)];
    mg = deg2rad(p.sun.marginDeg);       % back-off inside the band
    lbg = [lbg; cos(p.sun.thetaMax - mg)*ones(N*d, 1)];
    ubg = [ubg; cos(p.sun.thetaMin + mg)*ones(N*d, 1)];
end
blk.g = g; blk.lbg = lbg; blk.ubg = ubg;
blk.J = p.Tmax*sum(h.*sum(U, 1))/p.ocp.Jscale;
end
