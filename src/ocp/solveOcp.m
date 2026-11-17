function [sol, info] = solveOcp(ocp, w0, p, logFile, lam0)
%SOLVEOCP Run IPOPT from w0 (with the NLP parameters ocp.par, if any); IPOPT's
%   iteration log goes to logFile. lam0 (optional): struct with the initial
%   multipliers lam_x0, lam_g0 (used by a solver built with warm-start
%   options). sol.lamX, sol.lamG hold all multipliers (CasADi convention,
%   grad f + J_g' lam_g + lam_x = 0).
args = struct('x0', w0, 'lbx', ocp.lbx, 'ubx', ocp.ubx, ...
              'lbg', ocp.lbg, 'ubg', ocp.ubg);
if ~isempty(ocp.par), args.p = ocp.par; end
if nargin > 4 && ~isempty(lam0)
    args.lam_x0 = lam0.lam_x0; args.lam_g0 = lam0.lam_g0;
end
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
sol.lamG = full(res.lam_g);
sol.lamSun = sol.lamG(ocp.iSun);        % multipliers of (48), Problems 2-4
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
