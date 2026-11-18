function printProblem5(F, M, p)
%PRINTPROBLEM5 Console summary of Problem 5: the solve, Pi5, the
%   collision verification (i)-(iv) and one row per spacecraft.
fprintf('\n===== Problem 5 solution: Problem 4 + collision avoidance, dmin = %.1f m (NLP %.1f m) =====\n', ...
        p.collision.dmin, p.collision.dmin + p.collision.margin);
fprintf('Joint solver         : %s (%d iterations, %.1f s; start %d of %d; NLP built in %.1f s)\n', ...
        F.info.status, F.info.iter, F.info.wallTime, F.bestStart, numel(F.starts), F.buildTime);
fprintf('Total impulse J5*    : %.4f N s (J4* = %.4f N s): Pi5 = %.4f N s (%.3f%%)\n', ...
        M.J5, M.J4, M.Pi5, 100*M.Pi5rel);
fprintf('tau*, phiTau*        : %.1f s (%+.1f s), %.3f deg (%+.3f deg vs Problem 4)\n', ...
        M.tau, M.dTau, rad2deg(M.phiTau), M.dPhiDeg);
fprintf('dJ5/ddmin            : %.4f N s/m (multipliers), %.4f (central differences), %.4f (sweep slope)\n', ...
        M.dJddmin, M.dJddminFD, M.sweepSlopeAtDmin);
fprintf('Active pairs         : %s\n', mat2str(M.activePairs));
for a = M.activeIntervals
    fprintf('   pair (%d,%d) active on [%.0f, %.0f] s (%d collocation points)\n', a.pair, a.t0, a.t1, a.nPoints);
end
fprintf('Min. separation      : %.4f m (spacecraft %d and %d at t = %.1f s, propagated, exact)\n', ...
        M.minSep, M.minSepPair, M.minSepTime);
fprintf('\n i spacecraft            J_5,i    J_4,i   dm(kg)   arr.    abs  acc\n');
for i = 1:numel(F.res)
    r = F.res{i};
    fprintf('%2d %-18s %9.4f %8.4f %8.4f %6s %6s %4s\n', i, p.fleet.name(i), M.Ji(i), M.J4i(i), ...
            M.dmi(i), pf(r.term.passAll), pf(r.absPass), pf(r.accepted));
end
fprintf('\n dmin    J5(dmin)   dJ/ddmin  active  status\n');
for s = M.sweep
    fprintf('%5.1f %11.4f %10.4f %7d  %s\n', s.dmin, s.J_Ns, s.dJddmin, s.nActivePairs, s.status);
end
for f = fieldnames(M.checks)'
    fprintf('Check %-38s: %s\n', f{1}, pf(M.checks.(f{1})));
end
fprintf('(J tolerance of the cost checks: %.3f N s)\n', M.tolJ);
end

function s = pf(b)
if b, s = 'pass'; else, s = 'FAIL'; end
end
