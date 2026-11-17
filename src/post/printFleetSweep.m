function printFleetSweep(S, p)
%PRINTFLEETSWEEP Table of the Problem 4 sweep J_i(tf, phi) (kept class:
%   n/s for the north/south path, * where the -z_c guess won, r where a
%   retry was kept) and J4(phi),
%   and the per-spacecraft Problem 3 optima.
nS = size(S.Jbest, 2);
fprintf('\n===== Problem 4 sweep, tau = tf (J in N s; class n/s; * -z_c guess won, r retry) =====\n');
fprintf('%9s', 'phi(deg)');
for i = 1:nS, fprintf(' %15s', sprintf('J_%d', i)); end
fprintf(' %11s\n', 'J4');
cl = 'sn';                              % class -1 south, +1 north
for m = 1:numel(S.phis)
    tag = ''; if S.isPresc(m), tag = ' (P2)'; end
    fprintf('%9.1f', rad2deg(S.phis(m)));
    for i = 1:nS
        c = cl(max(S.clsBest(m, i), 0) + 1);
        star = ' '; if S.vBest(m, i) == 2, star = '*'; elseif S.vBest(m, i) > 2, star = 'r'; end
        fprintf(' %12.3f %c%c', S.Jbest(m, i), c, star);
    end
    mark = ''; if any(S.gridBest == m), mark = '  <- joint start'; end
    fprintf(' %11.3f%s%s\n', S.J4(m), tag, mark);
end
nFail = sum(~S.success(:));
fprintf('Failed subproblem solves: %d of %d guesses; %d (phi, i) retried, %d still failed\n', ...
        nFail, numel(S.success), numel(unique([[S.retry.sc]', [S.retry.phiDeg]'], 'rows')), ...
        sum(isnan(S.Jbest(:))));
fprintf('\n===== Problem 3 per spacecraft (J*_3,i) =====\n');
for i = 1:nS
    q = S.p3(i);
    fprintf('%d %-16s J*_3 %9.4f N s at phi %7.2f deg, tau %.1f s, class %+d, from %s (multi-start best %.4f)\n', ...
            i, p.fleet.name(i), q.J, rad2deg(q.phiTau), q.tau, q.cls, q.best, q.multistartBestJ);
end
fprintf('Sum J*_3,i = %.4f N s (lower bound of J4*); min_phi J4(phi) on the grid = %.4f N s\n', ...
        sum([S.p3.J]), min(S.J4(~S.isPresc)));
fprintf('Sweep wall time %.1f min, solver time %.1f min\n\n', S.wallTime/60, S.cpuTime/60);
end
