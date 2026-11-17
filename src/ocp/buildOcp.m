function ocp = buildOcp(p, ref, problem, parametric)
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
%   parametric = true: the initial state and the target point are NLP
%   parameters, ocp.par = [x_i,0./S; r_fp] from p, and setSpacecraft swaps
%   them for another spacecraft of the fleet (Problem 4: one compiled
%   solver per worker). Default false: they are constants, as in Problems
%   1-3. The two NLPs are equal up to rounding, but IPOPT's path on this
%   nonconvex problem is sensitive to it: the parametric form reaches
%   other nearby local optima of Problems 1 and 3 (measured: J1 168.80
%   instead of 166.37 N s), so Problems 1-3 keep the constant form.
%   The rows are built by ocpBlock, which Problem 4 (buildFleetOcp) shares.
import casadi.*
if nargin < 3, problem = 1; end
if nargin < 4, parametric = false; end
c = ocpFunctions(p, ref, problem >= 2);
N = c.N; d = c.d; nx = c.nx; nu = c.nu;

% Decision variables
nW = 1 + nx*(N + 1) + nx*N*d + nu*N + (problem >= 3);
w = MX.sym('w', nW);
parVal = [initialRelativeState(p)./c.S; p.rfp];
if parametric, par = MX.sym('par', nx + 3); else, par = parVal; end
ix.tau = 1;
ix.X   = 1 + (1:nx*(N + 1));
ix.XC  = ix.X(end) + (1:nx*N*d);
ix.U   = ix.XC(end) + (1:nu*N);
ix.phi = [];
if problem >= 3, ix.phi = nW; end

tau = p.tf*w(ix.tau);
if problem >= 3
    phiTau = pi*w(ix.phi);                    % decision variable, (56)
else
    phiTau = p.target.phase0 + p.target.rate*tau;   % prescribed, (43)
end
blk = ocpBlock(c, p, tau, phiTau, w(ix.X), w(ix.XC), w(ix.U), par);

% Bounds: (45) with the numerical lower bound, (44)
lbx = -inf(nW, 1); ubx = inf(nW, 1);
lbx(ix.tau) = p.tauMin/p.tf; ubx(ix.tau) = 1;
lbx(ix.U) = 0; ubx(ix.U) = 1;
lbx(ix.phi) = -1; ubx(ix.phi) = 1;          % (56), phiTau/pi

opts = struct('expand', true, 'print_time', false, 'ipopt', p.ocp.ipopt);
if problem >= 3, opts.ipopt.max_wall_time = p.ocp.p3MaxWallTime; end
ocp.nlp = struct('x', w, 'f', blk.J, 'g', blk.g);
ocp.par = [];
if parametric, ocp.nlp.p = par; ocp.par = parVal; end
ocp.opts = opts;
ocp.solver = nlpsol('solver', 'ipopt', ocp.nlp, opts);
ocp.physKey = c.physKey;
ocp.lbx = lbx; ocp.ubx = ubx;
ocp.lbg = blk.lbg; ocp.ubg = blk.ubg; ocp.problem = problem;
ocp.iSun = blk.iSun; ocp.iTerm = blk.iTerm;
ocp.ix = ix; ocp.N = N; ocp.d = d; ocp.nx = nx; ocp.nu = nu;
ocp.sigma = c.sigma; ocp.sNode = c.sNode; ocp.sColl = c.sColl; ocp.cRoot = c.cRoot;
ocp.C = c.C; ocp.D = c.D; ocp.S = c.S; ocp.F = c.F;
ocp.nW = nW; ocp.nG = numel(blk.g);
end
