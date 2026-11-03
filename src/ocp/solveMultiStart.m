function [sol, info, starts] = solveMultiStart(ocp, p)
%SOLVEMULTISTART Problem 3 from scratch (Sec. 5.2): one solve per start
%   arrival phase in p.ocp.p3Starts, each from the common initial guess
%   with the boundary states of that phase, and the best
%   Solve_Succeeded solution kept. starts lists every start.
phi0 = p.ocp.p3Starts;
starts = struct('phiStartDeg', {}, 'phiTauDeg', {}, 'J_Ns', {}, 'tau_s', {}, ...
                'status', {}, 'iterations', {}, 'time_s', {});
best = [];
for m = 1:numel(phi0)
    [s, in] = solveOcp(ocp, initialGuess(ocp, p, phi0(m)), p);
    starts(m) = struct('phiStartDeg', rad2deg(phi0(m)), 'phiTauDeg', rad2deg(s.phiTau), ...
        'J_Ns', s.J, 'tau_s', s.tau, 'status', in.status, 'iterations', in.iter, ...
        'time_s', in.wallTime);
    fprintf('Problem 3 start %d/%d: phiTau %7.1f -> %7.2f deg, J %.4f N s, %s (%d it., %.1f s)\n', ...
            m, numel(phi0), starts(m).phiStartDeg, starts(m).phiTauDeg, s.J, in.status, ...
            in.iter, in.wallTime);
    if in.success && (isempty(best) || s.J < best.sol.J)
        best = struct('sol', s, 'info', in, 'm', m);
    end
end
if isempty(best)
    error('LIFE:noStart', 'No Problem 3 start reached Solve_Succeeded.');
end
sol = best.sol; info = best.info;
info.bestStart = best.m;
info.multiStartTime = sum([starts.time_s]);
fprintf('Problem 3: best start %d, phiTau* = %.2f deg, J = %.4f N s\n', ...
        best.m, rad2deg(sol.phiTau), sol.J);
end
