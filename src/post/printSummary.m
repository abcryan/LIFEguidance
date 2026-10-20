function printSummary(res, p)
%PRINTSUMMARY Report of Sec. 4.3: arrival time, impulse, terminal errors.
t = res.term;
pf = @(b) char(ifelse(b));
fprintf('\n===== Final solution (level %d, N = %d) =====\n', res.level, res.sol.N);
fprintf('Solver status        : %s (%d iterations)\n', res.info.status, res.info.iter);
fprintf('Arrival time tau     : %.3f s (upper bound tf = %g s active: %d, multiplier %.3e)\n', ...
        res.sol.tau, p.tf, res.info.tauAtUpper, res.info.tauMultiplier);
fprintf('Total impulse        : %.4f N s  (dm = %.4f kg at Isp %g s)\n', res.J, res.dm, p.Isp);
fprintf('Max thrust           : %.4f N (Tmax %g N), bound violation %.2e N\n', ...
        res.maxThrust, p.Tmax, res.boundViolation);
fprintf('Terminal |r|         : %.3e m     (<= %g)  %s\n', t.r, p.gate.r, pf(t.pass.r));
fprintf('Terminal |v|         : %.3e m/s   (<= %g)  %s\n', t.v, p.gate.v, pf(t.pass.v));
fprintf('Terminal angle (60)  : %.3e rad   (<= %g)  %s   [atan2 form %.3e]\n', ...
        t.thetaAcos, p.gate.theta, pf(t.pass.theta), t.thetaAtan);
fprintf('Terminal |w|         : %.3e 1/s   (<= %g)  %s\n', t.w, p.gate.w, pf(t.pass.w));
fprintf('Doubling N 200->400  : dJ = %.3e (< %g)\n', res.doublingDJ, p.mesh.relTolJ);
fprintf('Last-level dJ        : %.3e\n', res.history.dJ(res.final));
fprintf('Prop. vs collocation : max |dr| %.2e m, |dv| %.2e m/s, dth %.2e rad, |dw| %.2e 1/s\n', ...
        res.disc.max.r, res.disc.max.v, res.disc.max.theta, res.disc.max.w);
fprintf('Absolute-model check : max |dr| %.2e m, |dv| %.2e m/s, dth %.2e rad, |dw| %.2e 1/s  %s\n', ...
        res.absDisc.max.r, res.absDisc.max.v, res.absDisc.max.theta, res.absDisc.max.w, pf(res.absPass));
fprintf('max ||q|-1|          : collocation %.2e, propagation %.2e; ||R''R-I|| %.2e\n', ...
        res.qNormColl, res.qNormProp, res.RorthProp);
fprintf('Burns                : %d, features converged across last two levels: %d\n', ...
        numel(res.burns), res.burnsConverged);
for b = 1:numel(res.burns)
    B = res.burns(b);
    fprintf('  burn %d: t = [%9.2f, %9.2f] s, peak %.3f N, impulse %.3f N s, |int f dt|_c %.3f N s, |int tau dt|_c %.3f N m s\n', ...
            b, B.t0, B.t1, B.peakT, B.impulse, norm(B.Fc), norm(B.Tc));
end
fprintf('ACCEPTED (status, (59)-(60), dJ < 1%%): %d\n\n', res.accepted);
end

function s = ifelse(b)
if b, s = "pass"; else, s = "FAIL"; end
end
