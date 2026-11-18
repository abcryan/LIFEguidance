function M = problem5Metrics(F, P4, p, feas)
%PROBLEM5METRICS Fleet quantities and the collision verification of
%   Problem 5 (Secs. 4.6, 5.2) from the validated solution F
%   (solveProblem5), the Problem 4 solution P4 (fleet, metrics) and the
%   formation check feas (checkCollisionFeasibility, step (i)):
%   - impulse per spacecraft, Pi5 = J5* - J4*, changes of tau* and phiTau*;
%   - (ii) collocation: every sampled separation >= dmin + margin, nonzero
%     multipliers only at the bound; the active pairs and intervals;
%   - (iii) propagation: every local minimum of d_ij(t) of all ten pairs,
%     also those not constrained in the NLP (p.collision.pairs =
%     "crossing"), located exactly (separationMinima), the smallest >= dmin;
%   - (iv) J4* <= J5*, J5(dmin) non-decreasing along the sweep and equal
%     to J4* below the Problem 4 minimum separation (and the solve at the
%     smallest distance from the Problem 4 optimum equal to J4*), dJ5/ddmin
%     from the multipliers against central differences;
%   - every spacecraft accepted with the unchanged gates.
%   CHOICE: "equal" and "non-decreasing" up to tolJ = 1e-4 J4* (0.14 N s),
%   which covers the re-solve scatter of up to 0.05 N s measured for the
%   Problem 4 subproblems; the derivative check uses the 10 % rule of
%   Problem 4.
cp = p.collision; nS = numel(F.res);
M.J5 = F.info.J; M.J4 = P4.fleet.info.J;
M.Ji = cellfun(@(r) r.J, F.res); M.J4i = cellfun(@(r) r.J, P4.fleet.res);
M.dmi = M.Ji/(p.g0*p.Isp); M.dm = sum(M.dmi);
M.Pi5 = M.J5 - M.J4; M.Pi5rel = M.Pi5/M.J4;
M.tau = F.info.tau; M.phiTau = F.info.phiTau;
M.dTau = M.tau - P4.fleet.info.tau;
M.dPhiDeg = rad2deg(angle(exp(1i*(M.phiTau - P4.fleet.info.phiTau))));
M.cls = cellfun(@(r) pathClass(r.sun.azimuth), F.res);
M.accepted = cellfun(@(r) r.accepted, F.res);
M.feasibility = feas;
tolJ = 1e-4*M.J4; M.tolJ = tolJ;

% (ii) collocation solution
cs = F.colloc;
M.colloc = struct('minPerPair', cs.minPerPair, 'bound', cs.bound, ...
                  'slackAtActiveM', cs.slackAtActive);
act = cs.active; pairs = F.pairs;
M.activePairs = pairs(any(act, 2), :);
M.activeIntervals = struct('pair', {}, 't0', {}, 't1', {}, 'nPoints', {});
for k = find(any(act, 2))'
    on = act(k, :); de = diff([0 on 0]); a = find(de == 1); b = find(de == -1) - 1;
    for m = 1:numel(a)
        M.activeIntervals(end + 1) = struct('pair', pairs(k, :), 't0', cs.t(a(m)), ...
            't1', cs.t(b(m)), 'nPoints', b(m) - a(m) + 1);
    end
end

% (iii) propagated trajectories: every local minimum of every pair (also
% the pairs not constrained in the NLP), exactly
allPairs = nchoosek(1:nS, 2);
M.minima = separationMinima(F.res, p, allPairs);
[M.minSep, km] = min([M.minima.dMin]);
M.minSepPair = M.minima(km).pair; M.minSepTime = M.minima(km).tMin;

% Distances of all pairs on the common propagation grid (figures), and
% Problem 4's
t = F.res{1}.prop.t;
D = zeros(size(allPairs, 1), numel(t));
for k = 1:size(allPairs, 1)
    D(k, :) = vecnorm(F.res{allPairs(k,1)}.prop.xP(8:10,:) - F.res{allPairs(k,2)}.prop.xP(8:10,:));
end
M.constrainedPairs = pairs;
M.dist = struct('t', t, 'pairs', allPairs, 'D', D, 'min', M.minSep, 'minPair', M.minSepPair, ...
                'minTime', M.minSepTime);
M.dist4 = P4.metrics.dist;

% (iv) costs and sensitivity
sw = F.sweep; ok = [sw.success];
M.sweep = sw;
below = [sw.dmin] + cp.margin < P4.metrics.dist.min;
M.dJddmin = F.sens.dJddmin; M.dJddminFD = F.fd.dJddmin;
c.formations_satisfy_dmin = feas.ok;
c.collocation_separations_ge_bound = min(cs.minPerPair) >= cs.bound - 1e-3;
c.multipliers_only_at_bound = cs.slackAtActive <= 1e-3;
c.propagated_min_separation_ge_dmin = M.minSep >= cp.dmin;
c.J4_le_J5 = M.J4 <= M.J5 + tolJ;
c.sweep_converged = all(ok);
Jok = [sw(ok).J_Ns];
c.sweep_nondecreasing = all(diff(Jok) >= -tolJ);
c.J5_equals_J4_below_P4_min_separation = any(below & ok) && ...
    all(abs([sw(below & ok).J_Ns] - M.J4) <= tolJ);
c.J5_low_from_P4_equals_J4 = F.lowFromP4.success && abs(F.lowFromP4.J_Ns - M.J4) <= tolJ;
c.sensitivity_matches_finite_differences = all(F.fd.success) && ...
    abs(M.dJddmin - M.dJddminFD) <= 0.1*max(abs([M.dJddmin, M.dJddminFD]));
c.solve_succeeded = F.info.success;
c.all_spacecraft_accepted = all(M.accepted);
M.checks = c;
% Expectations of Sec. 5.1 (reported, not enforced)
crossPairs = reshape([F.crossing.pair], 2, [])';
e.only_crossing_pairs_active = isequal(sortrows(M.activePairs), sortrows(crossPairs));
M.expectations = e;
M.sweepSlopeAtDmin = NaN;
swOk = sw(ok); i0 = find([swOk.dmin] == cp.dmin);      % converged neighbours only
if ~isempty(i0) && i0 > 1 && i0 < numel(swOk)
    M.sweepSlopeAtDmin = (swOk(i0+1).J_Ns - swOk(i0-1).J_Ns)/(swOk(i0+1).dmin - swOk(i0-1).dmin);
end
M.nS = nS;
end
