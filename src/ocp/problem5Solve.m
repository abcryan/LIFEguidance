function r = problem5Solve(fo, p, w0, dmins, sides, pairs)
%PROBLEM5SOLVE One Problem 5 solve, or a chain of them, on the Problem 5
%   NLP fo (buildFleetOcp with collision rows), in the client or on a pool
%   worker (problem5Worker, with the worker's own compiled NLP):
%   - w0: primal start (the Problem 4 layout);
%   - dmins: one collision distance, or a continuation chain, each solve
%     started from the previous one (from w0 again after a failed step);
%   - sides, pairs (optional): before the first solve, displace the two
%     collectors of each pairs(m,:) in opposite directions along
%     sides(m)*z_* by (dmin/2) sin(pi nu) (displaceStart).
%   r(k) per distance: dmin, the single-spacecraft solutions and the joint
%   info (with x and lamG), the collision sensitivity and the sampled
%   separations (collisionSensitivity, collocationSeparations).
if nargin >= 6 && ~isempty(pairs)
    w0 = displaceStart(fo, w0, pairs, sides, dmins(1), p);
end
mg = p.collision.margin;
w = w0;
for k = numel(dmins):-1:1, r(k).dmin = dmins(k); end     % preallocate
for k = 1:numel(dmins)
    fk = setCollision(fo, dmins(k), mg);
    [sols, in] = solveFleetOcp(fk, fo.solver, w, [], p);
    r(k).sols = sols; r(k).info = in;
    r(k).sens = collisionSensitivity(in, fo, dmins(k), mg, p);
    r(k).cs = collocationSeparations(sols, fo, in, dmins(k), mg);
    if in.success, w = in.x; end
end
end

function w = displaceStart(fo, w, pairs, sides, dmin, p)
% Collectors of each crossing pair displaced in opposite directions along
% +-z_* by (dmin/2) sin(pi nu), with the matching velocity, at the nodes
% and collocation points (scaled states).
tau = p.tf*w(1);
S = fo.S;
v1 = fo.view(1);
nuN = v1.sNode(:)'; nuC = v1.sColl(:)';
for m = 1:size(pairs, 1)
    for e = 1:2
        i = pairs(m, e); sgn = sides(m)*(3 - 2*e);    % + for the first, - for the second
        v = fo.view(i);
        X = reshape(w(v.ix.X), fo.nx, []); XC = reshape(w(v.ix.XC), fo.nx, []);
        X(10, :) = X(10, :) + sgn*dmin/2*sin(pi*nuN)/S(10);
        X(13, :) = X(13, :) + sgn*dmin/2*pi/tau*cos(pi*nuN)/S(13);
        XC(10, :) = XC(10, :) + sgn*dmin/2*sin(pi*nuC)/S(10);
        XC(13, :) = XC(13, :) + sgn*dmin/2*pi/tau*cos(pi*nuC)/S(13);
        w(v.ix.X) = X(:); w(v.ix.XC) = XC(:);
    end
end
end

function sn = collisionSensitivity(in, fo, dk, mg, p)
% dJ/ddmin = sum_k dJ/dlbg_k * dlbg_k/ddmin with dJ/dlbg = -lam_g (CasADi
% sign convention: lam_g <= 0 on an active lower bound) and
% lbg = ((dmin + margin)/scale)^2; J scaled by Jscale in the NLP.
lam = in.lamG(fo.iCol);
dlbg = 2*(dk + mg)/fo.colScale^2;
sn.dJddmin = p.ocp.Jscale*sum(-lam)*dlbg;
act = abs(lam) > activeThreshold(lam);
sn.activePairs = unique(fo.colPair(act))';
sn.lamMax = max(abs(lam));
end

function cs = collocationSeparations(sols, fo, in, dk, mg)
% Separations of every pair at the collocation points (SI); the active
% rows (multiplier above threshold) and the largest distance of an active
% row from the imposed bound (complementarity).
nP = size(fo.pairs, 1); nC = fo.N*fo.d;
D = zeros(nP, nC);
for k = 1:nP
    D(k, :) = vecnorm(sols{fo.pairs(k,1)}.XC(8:10, :) - sols{fo.pairs(k,2)}.XC(8:10, :));
end
cs.minPerPair = min(D, [], 2)';
cs.D = D;
cs.t = sols{1}.tColl;
lam = reshape(in.lamG(fo.iCol), nC, nP)';
cs.bound = dk + mg;
cs.active = abs(lam) > activeThreshold(lam);
cs.lam = lam;
cs.slackAtActive = max([0; abs(D(cs.active) - cs.bound)]);   % [m]
end

function thr = activeThreshold(lam)
% CHOICE: a collision row is active if its multiplier exceeds 1e-4 of the
% largest. Measured (N = 10 check): the rows at the bound carry 3e-2-1e-1,
% all others <= 3e-6 at >= 3 mm from it (IPOPT's barrier residue, lambda
% times slack ~ const), four orders of magnitude apart. At N = 400 the
% multipliers fall off smoothly at the ends of an active arc instead.
% Absolute floor 1e-6: when no row is active (e.g. dmin below the Problem 4
% separation), every multiplier is barrier residue and none counts.
thr = max(1e-4*max(abs(lam(:))), 1e-6);
end
