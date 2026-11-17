function M = fleetMetrics(F, S, p)
%FLEETMETRICS Fleet quantities and checks of Problem 4 (Secs. 4.5, 4.6,
%   5.2, 5.3) from the validated joint solution F (solveFleet) and the
%   sweep S (fleetSweep.m):
%   - impulse and propellant per spacecraft, the combiner's share;
%   - the price of the shared variables Pi4 = J4* - sum_i J*_3,i, eq. (66);
%   - consistency checks (local optima, so checks only): J*_3,i <= the
%     minimum of spacecraft i's own sweep, sum_i J*_3,i <= J4*,
%     J_4,i >= J*_3,i, J4* <= the sweep minimum (grid), J4* <= sum_i J_i at
%     Problem 2's prescribed phase, sensitivities summing to zero at an
%     interior phase and agreeing with the finite differences;
%   - the attitude-path class of every spacecraft;
%   - the distances between all pairs over time (reported, not constrained);
%   - the free-space expectations of Sec. 5.1, reported against the result.
nS = numel(F.res);
r1 = F.res{1};
M.J4 = F.info.J;
M.Ji = cellfun(@(r) r.J, F.res);
M.dmi = M.Ji/(p.g0*p.Isp);
M.dm = sum(M.dmi);
M.share = M.Ji/M.J4;
M.tau = r1.sol.tau; M.phiTau = r1.sol.phiTau;
M.J3i = [S.p3.J]; M.phi3i = [S.p3.phiTau];
M.Pi4 = M.J4 - sum(M.J3i);
M.Pi4rel = M.Pi4/M.J4;
J4g = S.J4; J4g(S.isPresc) = NaN;
[M.sweepMin, k] = min(J4g); M.sweepMinPhi = S.phis(k);
M.prescBound = S.J4(S.isPresc);
M.dJdphi = F.info.dJdphi;
M.dJdphiFD = F.fd.dJdphi;
M.cls = cellfun(@(r) pathClass(r.sun.azimuth), F.res);
M.accepted = cellfun(@(r) r.accepted, F.res);

tol = 1e-9;
c.sum_J3_le_J4 = sum(M.J3i) <= M.J4*(1 + tol);
% Problem 3 per spacecraft at least as good as its own sweep (Sec. 5.2)
M.sweepMinI = min(S.Jbest(~S.isPresc, :), [], 1);
c.each_J3i_le_own_sweep_min = all(M.J3i <= M.sweepMinI*(1 + 1e-6));
c.each_J4i_ge_J3i = all(M.Ji >= M.J3i*(1 - tol));
c.J4_le_sweep_min = M.J4 <= M.sweepMin*(1 + 1e-6);
c.J4_le_prescribed_phase_bound = M.J4 <= M.prescBound*(1 + 1e-6);
M.sensSum = sum(M.dJdphi);
c.phase_interior = ~F.info.phiAtBound;
c.sensitivities_sum_to_zero = abs(M.sensSum) <= 1e-6*max(1, sum(abs(M.dJdphi)));
M.fdErr = abs(M.dJdphiFD - M.dJdphi);
c.finite_differences_agree = all(F.fd.success) && max(M.fdErr) <= 0.1*max(abs(M.dJdphi));
c.solve_succeeded = F.info.success;
c.all_spacecraft_accepted = all(M.accepted);
M.checks = c;

% Distances between all pairs on the common propagation grid
t = r1.prop.t;
for i = 2:nS
    assert(isequal(F.res{i}.prop.t, t), 'Propagation grids differ.');
end
pairs = nchoosek(1:nS, 2);
D = zeros(size(pairs, 1), numel(t));
for k = 1:size(pairs, 1)
    D(k,:) = vecnorm(F.res{pairs(k,1)}.prop.xP(8:10,:) - F.res{pairs(k,2)}.prop.xP(8:10,:));
end
[dmin, kk] = min(D(:));
[kp, kt] = ind2sub(size(D), kk);
M.dist = struct('t', t, 'pairs', pairs, 'D', D, 'min', dmin, 'minPair', pairs(kp,:), ...
                'minTime', t(kt), 'atStart', D(:,1)', 'atArrival', D(:,end)');

% Free-space expectations of Sec. 5.1 (reported, not enforced)
e.tau_at_tf = F.info.tauAtUpper;
e.phi_within_15_deg_of_0 = abs(rad2deg(M.phiTau)) <= 15;
e.Pi4_below_1_percent = M.Pi4rel < 0.01;
e.combiner_share_above_half = M.share(5) > 0.5;
% individually preferred phases of the diagonal pairs (1,4) and (2,3)
pm = @(a) angle(mean(exp(1i*a)));
e.diagonal_pair_phase_difference_deg = rad2deg(angle(exp(1i*(pm(M.phi3i([1 4])) - pm(M.phi3i([2 3]))))));
M.expectations = e;
end
