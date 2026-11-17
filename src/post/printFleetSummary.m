function printFleetSummary(F, M, p)
%PRINTFLEETSUMMARY Report of Problem 4: one row per spacecraft, the fleet
%   totals, the price of the shared variables, the checks and the
%   free-space expectations.
nS = numel(F.res);
pf = @(b) char(string(ifelse(b)));
in = F.info;
fprintf('\n===== Problem 4 solution: %d spacecraft, shared tau and phiTau (N = %d) =====\n', nS, p.ocp.N);
fprintf('Joint solver         : %s (%d iterations, %.1f s; best of %d starts, %.1f s in all; NLP built in %.1f s)\n', ...
        in.status, in.iter, in.wallTime, numel(F.starts), F.jointTime, F.buildTime);
fprintf('Arrival time tau*    : %.3f s (tf = %g s active: %d, multiplier %.3e)\n', M.tau, p.tf, ...
        in.tauAtUpper, in.tauMultiplier);
fprintf('Arrival phase phi*   : %.3f deg (free, shared; phase multiplier %.1e)\n', rad2deg(M.phiTau), ...
        in.phiMultiplier);
fprintf('Total impulse J4*    : %.4f N s (dm = %.4f kg at Isp %g s)\n', M.J4, M.dm, p.Isp);
fprintf('Sum of J*_3,i        : %.4f N s -> price of the shared variables Pi4 = %.4f N s (%.3f%% of J4*)\n', ...
        sum(M.J3i), M.Pi4, 100*M.Pi4rel);
fprintf('Sweep minimum J4(phi): %.4f N s at %.1f deg; Problem 2 phase bound %.4f N s\n', ...
        M.sweepMin, rad2deg(M.sweepMinPhi), M.prescBound);
fprintf('\n%2s %-16s %10s %6s %8s %10s %8s %9s %9s %5s %8s %6s %6s %4s\n', 'i', 'spacecraft', ...
        'J_4,i', 'share', 'dm(kg)', 'J*_3,i', 'phi*_3,i', 'dJ/dphi', 'FD', 'class', 'Sun(deg)', 'arr.', 'abs', 'acc');
cl = ["south", "", "north"];
for i = 1:nS
    r = F.res{i};
    fprintf('%2d %-16s %10.4f %5.1f%% %8.4f %10.4f %8.2f %9.3f %9.3f %5s %8.4f %6s %6s %4s\n', ...
            i, p.fleet.name(i), M.Ji(i), 100*M.share(i), M.dmi(i), M.J3i(i), rad2deg(M.phi3i(i)), ...
            M.dJdphi(i), M.dJdphiFD(i), cl(M.cls(i) + 2), r.sun.maxViolDeg, pf(r.term.passAll), ...
            pf(r.absPass), pf(r.accepted));
end
fprintf(['(dJ/dphi in N s/rad from the terminal multipliers, FD: slope of a quadratic fit through ' ...
         're-solves over +-%g deg; Sun: max violation)\n'], 2*F.fd.stepDeg);
fprintf('Sum of sensitivities : %.3e N s/rad\n', M.sensSum);
for i = 1:nS
    t = F.res{i}.term;
    fprintf('  sc %d arrival: |dr| %.2e m, |dv| %.2e m/s, xi %.2e rad, |dw| %.2e 1/s; Sun active at %d points\n', ...
            i, t.r, t.v, t.thetaAcos, t.w, numel(F.res{i}.sun.active));
end
d = M.dist;
fprintf('Min. distance        : %.1f m between spacecraft %d and %d at t = %.0f s (not constrained)\n', ...
        d.min, d.minPair, d.minTime);
fn = fieldnames(M.checks);
for k = 1:numel(fn)
    fprintf('Check %-30s: %s\n', fn{k}, pf(M.checks.(fn{k})));
end
e = M.expectations;
fprintf('Expectations (Sec. 5.1, free-space estimate): tau* = tf %d, |phi*| <= 15 deg %d, Pi4 < 1%% %d, combiner > 50%% %d, diagonal pairs %.1f deg apart (doc ~36)\n\n', ...
        e.tau_at_tf, e.phi_within_15_deg_of_0, e.Pi4_below_1_percent, e.combiner_share_above_half, ...
        e.diagonal_pair_phase_difference_deg);
end

function s = ifelse(b)
if b, s = "pass"; else, s = "FAIL"; end
end
