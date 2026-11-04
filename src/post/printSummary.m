function printSummary(res, p)
%PRINTSUMMARY Report of Sec. 4.3: arrival time, impulse, terminal errors.
t = res.term;
pf = @(b) char(ifelse(b));
fprintf('\n===== Problem %d solution (N = %d uniform intervals) =====\n', res.problem, res.sol.N);
fprintf('Solver status        : %s (%d iterations, %.1f s)\n', res.info.status, ...
        res.info.iter, res.info.wallTime);
fprintf('Arrival time tau     : %.3f s (upper bound tf = %g s active: %d, multiplier %.3e)\n', ...
        res.sol.tau, p.tf, res.info.tauAtUpper, res.info.tauMultiplier);
fprintf('Total impulse        : %.4f N s  (dm = %.4f kg at Isp %g s)\n', res.J, res.dm, p.Isp);
if res.problem >= 3
    bnd = '';
    if res.info.phiAtBound, bnd = '; at its bound'; end
    fprintf('Arrival phase phiTau : %.3f deg (free%s), formation phase at t = 0: %.3f deg\n', ...
            rad2deg(res.sol.phiTau), bnd, rad2deg(res.sol.phase0));
else
    fprintf('Arrival phase phiTau : %.3f deg (prescribed, phi(0) = %.3f deg)\n', ...
            rad2deg(res.sol.phiTau), rad2deg(res.sol.phase0));
end
fprintf('Max thrust           : %.4f N (Tmax %g N), bound violation %.2e N\n', ...
        res.maxThrust, p.Tmax, res.boundViolation);
fprintf('Terminal |r|         : %.3e m     (<= %g)  %s\n', t.r, p.gate.r, pf(t.pass.r));
fprintf('Terminal |v|         : %.3e m/s   (<= %g)  %s\n', t.v, p.gate.v, pf(t.pass.v));
fprintf('Terminal angle xi    : %.3e rad   (<= %g)  %s   [atan2 form %.3e]\n', ...
        t.thetaAcos, p.gate.theta, pf(t.pass.theta), t.thetaAtan);
fprintf('Terminal |w|         : %.3e 1/s   (<= %g)  %s\n', t.w, p.gate.w, pf(t.pass.w));
fprintf('Prop. vs collocation : max |dr| %.2e m, |dv| %.2e m/s, dth %.2e rad, |dw| %.2e 1/s\n', ...
        res.disc.max.r, res.disc.max.v, res.disc.max.theta, res.disc.max.w);
fprintf('Absolute-model check : max |dr| %.2e m, |dv| %.2e m/s, dth %.2e rad, |dw| %.2e 1/s  %s\n', ...
        res.absDisc.max.r, res.absDisc.max.v, res.absDisc.max.theta, res.absDisc.max.w, pf(res.absPass));
fprintf('max ||q|-1|          : collocation %.2e, propagation %.2e; ||R''R-I|| %.2e\n', ...
        res.qNormColl, res.qNormProp, res.RorthProp);
sn = res.sun;
fprintf('Sun angle (48)       : anti-Sun angle of n^i in [%.2f, %.2f] deg (band [%.0f, %.0f]), max violation %.3f deg (<= %g)  %s\n', ...
        min(sn.antiDeg), max(sn.antiDeg), 180 - rad2deg(p.sun.thetaMax), 180 - rad2deg(p.sun.thetaMin), ...
        sn.maxViolDeg, p.gate.sunDeg, sunStatus(sn.pass, res.problem));
fprintf('                       time outside the band (> gate) %.0f s, at collocation points %.2e deg, formation-centre error %.1e rad\n', ...
        sn.timeOutside, sn.collMaxViolDeg, sn.centreErr);
for b = 1:numel(sn.intervals)
    I = sn.intervals(b);
    fprintf('  outside: t = [%7.0f, %7.0f] s, anti-Sun angle [%.1f, %.1f] deg\n', I.t0, I.t1, I.antiDegMin, I.antiDegMax);
end
if res.problem >= 2
    fprintf('Sun constraint active: %d collocation points (%s)\n', numel(sn.active), sn.activeBound);
    for a = sn.active
        fprintf('  touch: t = %8.1f s, anti-Sun angle %.4f deg (bound %.0f), multiplier %+.2e\n', ...
                a.t, a.antiDeg, a.bound, a.lambda);
    end
end
fprintf('Burns                : %d\n', numel(res.burns));
for b = 1:numel(res.burns)
    B = res.burns(b);
    fprintf('  burn %d: t = [%9.2f, %9.2f] s, peak %.3f N, impulse %.3f N s, |int f dt|_c %.3f N s, |int tau dt|_c %.3f N m s\n', ...
            b, B.t0, B.t1, B.peakT, B.impulse, norm(B.Fc), norm(B.Tc));
end
if res.problem >= 2
    fprintf('ACCEPTED (Solve_Succeeded, arrival checks, Sun angle): %d\n\n', res.accepted);
else
    fprintf('ACCEPTED (Solve_Succeeded and arrival checks): %d\n\n', res.accepted);
end
end

function s = sunStatus(pass, problem)
%SUNSTATUS Gate result for Problem 2; Problem 1 does not impose (48).
if problem >= 2, s = char(ifelse(pass));
elseif pass, s = 'within (not imposed)';
else, s = 'violated (not imposed)';
end
end

function s = ifelse(b)
if b, s = "pass"; else, s = "FAIL"; end
end
