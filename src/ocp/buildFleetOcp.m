function ocp = buildFleetOcp(p, ref, ipoptExtra, collision)
%BUILDFLEETOCP Problem 4, eqs. (57)-(64): the five spacecraft of the fleet
%   on a common normalized mesh, coupled only through the shared arrival
%   time tau and arrival phase phiTau. Each block is the Problem 3 block of
%   ocpBlock (collocation, initial state, terminal condition (39) with its
%   own target point, Sun-angle constraint) for spacecraftParams(p, i).
%   Decision vector w = [tau/tf; phiTau/pi; w_1; ...; w_nS] with
%   w_i = [X_i(:); XC_i(:); U_i(:)] as in buildOcp; the constraint Jacobian
%   is block-angular (two coupling columns). Objective: total fleet
%   impulse (57). ocp.view(i) is a single-spacecraft OCP struct indexing
%   into the full w (for packW, unpackW). ocp.dgTdphi(w) gives the
%   derivative of the terminal rows by phiTau/pi (sensitivities).
%   ipoptExtra: IPOPT options added to p.ocp.ipopt (e.g. warm start).
%   collision = true (Problem 5, eq. (67)): appends, for all ten pairs (or
%   only the pairs listed, collision = k x 2 matrix) and every
%   collocation point, the rows g = |r_i - r_j|^2/dmin^2 of the
%   pointing-frame positions (p.collision.dmin as fixed scale). Their lower
%   bound is set per solve (setCollision), so one solver serves every
%   dmin; lbg = 0 is Problem 4. No new variables: w is the Problem 4 layout.
%   ocp.iCol (rows), ocp.colPair (pair index per row), ocp.colTime
%   (normalized time per row), ocp.pairs (the constrained pairs, k x 2).
%   Without it the NLP is the
%   Problem 4 NLP, unchanged.
import casadi.*
if nargin < 3 || isempty(ipoptExtra), ipoptExtra = struct(); end
if nargin < 4, collision = false; end
nS = numel(p.fleet.point);
pk = arrayfun(@(i) spacecraftParams(p, i), 1:nS, 'UniformOutput', false);
c = ocpFunctions(p, ref, true);
N = c.N; d = c.d; nx = c.nx; nu = c.nu;
nB = nx*(N + 1) + nx*N*d + nu*N;
nW = 2 + nS*nB;
w = MX.sym('w', nW);
tau = p.tf*w(1);
phiTau = pi*w(2);

g = cell(nS, 1); lbg = cell(nS, 1); ubg = cell(nS, 1);
J = 0; nG = 0; iTerm = [];
lbx = -inf(nW, 1); ubx = inf(nW, 1);
for i = 1:nS
    ci = c;                              % spacecraft-specific physics
    if ~isequal(physicsKey(pk{i}), c.physKey), ci = ocpFunctions(pk{i}, ref, true); end
    off = 2 + (i - 1)*nB;
    v.ix.tau = 1; v.ix.phi = 2;
    v.ix.X  = off + (1:nx*(N + 1));
    v.ix.XC = v.ix.X(end) + (1:nx*N*d);
    v.ix.U  = v.ix.XC(end) + (1:nu*N);
    par = [initialRelativeState(pk{i})./c.S; pk{i}.rfp];
    blk = ocpBlock(ci, pk{i}, tau, phiTau, w(v.ix.X), w(v.ix.XC), w(v.ix.U), par);
    v.iSun = nG + blk.iSun; v.iTerm = nG + blk.iTerm;
    v.iG = nG + (1:numel(blk.g))';
    v.N = N; v.d = d; v.nx = nx; v.nu = nu; v.S = c.S; v.nW = nW;
    v.sigma = c.sigma; v.sNode = c.sNode; v.sColl = c.sColl; v.cRoot = c.cRoot;
    v.par = []; v.problem = 4;
    view(i) = v; %#ok<AGROW>
    g{i} = blk.g; lbg{i} = blk.lbg; ubg{i} = blk.ubg;
    J = J + blk.J;
    nG = nG + numel(blk.g);
    iTerm = [iTerm; v.iTerm]; %#ok<AGROW>
    lbx(v.ix.U) = 0; ubx(v.ix.U) = 1;    % (61)
end
g = vertcat(g{:});
lbg = vertcat(lbg{:}); ubg = vertcat(ubg{:});
% Collision avoidance (67): squared separations of the position states
% (rows 8:10 of XC, scaled by S(8:10)) at every collocation point
pairs = nchoosek(1:nS, 2);
if ~islogical(collision), pairs = collision; collision = ~isempty(pairs); end
ocp.pairs = pairs; ocp.iCol = []; ocp.colPair = []; ocp.colTime = [];
if collision
    nP = size(pairs, 1); nC = N*d;
    sc = c.S(8)/p.collision.dmin;        % positions in units of dmin
    pos = cell(1, nS);
    for i = 1:nS
        XCi = reshape(w(view(i).ix.XC), nx, nC);
        pos{i} = XCi(8:10, :);
    end
    gc = cell(nP, 1);
    for k = 1:nP
        dr = pos{pairs(k,1)} - pos{pairs(k,2)};
        gc{k} = sc^2*sum1(dr.^2)';
    end
    ocp.iCol = nG + (1:nP*nC)';
    ocp.colPair = kron((1:nP)', ones(nC, 1));
    ocp.colTime = repmat(c.sColl(:), nP, 1);
    g = [g; vertcat(gc{:})];
    lbg = [lbg; zeros(nP*nC, 1)]; ubg = [ubg; inf(nP*nC, 1)];
    nG = nG + nP*nC;
end
% (63) with the numerical lower bound, (64)
lbx(1) = p.tauMin/p.tf; ubx(1) = 1;
lbx(2) = -1; ubx(2) = 1;

opts = struct('expand', true, 'print_time', false, 'ipopt', p.ocp.ipopt);
opts.ipopt.max_wall_time = p.ocp.p4MaxWallTime;
if collision, opts.ipopt.max_wall_time = p.ocp.p5MaxWallTime; end
for f = fieldnames(ipoptExtra)'
    opts.ipopt.(f{1}) = ipoptExtra.(f{1});
end
ocp.nlp = struct('x', w, 'f', J, 'g', g);
ocp.opts = opts;
ocp.solver = nlpsol('fleet', 'ipopt', ocp.nlp, opts);
JT = jacobian(g(iTerm), w);             % w(2) is not a symbol by itself
ocp.dgTdphi = Function('dgTdphi', {w}, {JT(:, 2)});
ocp.iTermAll = iTerm;
ocp.view = view;
ocp.pk = pk;
ocp.nS = nS; ocp.nB = nB; ocp.nW = nW; ocp.nG = nG;
ocp.N = N; ocp.d = d; ocp.nx = nx; ocp.nu = nu; ocp.S = c.S;
ocp.lbx = lbx; ocp.ubx = ubx;
ocp.lbg = lbg; ocp.ubg = ubg;
ocp.collision = collision; ocp.colScale = p.collision.dmin;
end
