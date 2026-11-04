function [sol, info] = solveOcp(ocp, w0, p, logFile)
%SOLVEOCP Run IPOPT from w0; IPOPT's iteration log goes to logFile.
args = struct('x0', w0, 'lbx', ocp.lbx, 'ubx', ocp.ubx, ...
              'lbg', ocp.lbg, 'ubg', ocp.ubg);
if nargin > 3 && ~isempty(logFile)
    if isfile(logFile), delete(logFile); end
    diary(logFile); cleanup = onCleanup(@() diary('off'));
end
tic;
res = ocp.solver.call(args);
info.wallTime = toc;
st = ocp.solver.stats();
info.status = st.return_status;
info.success = strcmp(st.return_status, 'Solve_Succeeded');
info.iter = st.iter_count;
sol = unpackW(ocp, res.x, p);
sol.lamX = full(res.lam_x);
lamG = full(res.lam_g);
sol.lamSun = lamG(ocp.iSun);            % multipliers of (48), Problems 2-3
g = full(res.g);
eq = ocp.lbg == ocp.ubg;                % collocation, continuity, boundary
info.maxDefect = max(abs(g(eq) - ocp.lbg(eq)));
info.objScaled = full(res.f);
% Multiplier of the arrival-time bound (positive: upper bound active)
info.tauMultiplier = sol.lamX(ocp.ix.tau);
info.tauAtUpper = abs(sol.tau - p.tf) < 1e-6*p.tf;
info.boundViolation = max([0; -sol.U(:); sol.U(:) - 1])*p.Tmax;
% Problem 3: arrival phase at its bound -pi or pi (55)
info.phiAtBound = ~isempty(ocp.ix.phi) && abs(abs(sol.phiTau) - pi) < 1e-6;
end
