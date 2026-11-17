function [sols, info] = solveFleetOcp(fo, solver, w0, lam0, p)
%SOLVEFLEETOCP One solve of the joint Problem 4 NLP (buildFleetOcp) with
%   solver (fo.solver or one with other options) from w0 and, if lam0 is
%   not empty, the multipliers lam0.lam_x0, lam0.lam_g0. sols{i} is the
%   solution of spacecraft i in the single-spacecraft layout (unpackW of
%   view i, shared tau and phiTau), with its Sun-constraint multipliers.
%   info: status, iterations, time, total impulse, and the sensitivities
%   dJ_i/dphiTau [N s/rad] from the multipliers of the terminal rows (39)
%   of each block (envelope theorem; they sum to zero at an interior
%   phase optimum), with lam_x and lam_g for later warm starts.
args = struct('x0', w0, 'lbx', fo.lbx, 'ubx', fo.ubx, 'lbg', fo.lbg, 'ubg', fo.ubg);
if ~isempty(lam0)
    args.lam_x0 = lam0.lam_x0; args.lam_g0 = lam0.lam_g0;
end
tic;
res = solver.call(args);
info.wallTime = toc;
st = solver.stats();
info.status = st.return_status;
info.success = strcmp(st.return_status, 'Solve_Succeeded');
info.iter = st.iter_count;
x = full(res.x); lamX = full(res.lam_x); lamG = full(res.lam_g);
g = full(res.g);
eq = fo.lbg == fo.ubg;
info.maxDefect = max(abs(g(eq) - fo.lbg(eq)));
info.tauMultiplier = lamX(1);
info.phiMultiplier = lamX(2);
dg = full(fo.dgTdphi(x));
sols = cell(1, fo.nS);
info.dJdphi = zeros(1, fo.nS);
for i = 1:fo.nS
    v = fo.view(i);
    s = unpackW(v, x, fo.pk{i});
    s.lamSun = lamG(v.iSun);
    sols{i} = s;
    rows = ismember(fo.iTermAll, v.iTerm);
    info.dJdphi(i) = p.ocp.Jscale/pi*(lamG(v.iTerm)'*dg(rows));
end
info.J = sum(cellfun(@(s) s.J, sols));
info.Ji = cellfun(@(s) s.J, sols);
info.tau = sols{1}.tau; info.phiTau = sols{1}.phiTau;
info.tauAtUpper = abs(info.tau - p.tf) < 1e-6*p.tf;
info.phiAtBound = abs(abs(info.phiTau) - pi) < 1e-6;
Uall = cell2mat(cellfun(@(s) s.U(:), sols, 'UniformOutput', false)');
info.boundViolation = max([0; -Uall; Uall - 1])*p.Tmax;
info.x = x; info.lamX = lamX; info.lamG = lamG;
end
