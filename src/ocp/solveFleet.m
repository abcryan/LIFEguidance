function F = solveFleet(p, ref, S)
%SOLVEFLEET Problem 4 from the sweep S of fleetSweep.m (Sec. 5.2):
%   1. the joint NLP (buildFleetOcp, tau and phiTau free) from each of the
%      p.ocp.p4.nJointStarts best grid points of the sweep: the assembled
%      subproblem solutions as primal start, with continuationOptions so
%      that the joint solve keeps their local optima while it moves the
%      shared phase. Best Solve_Succeeded kept.
%   2. the sensitivities dJ_i/dphiTau from the terminal multipliers, and a
%      finite-difference check: each subproblem at phi* +- p.ocp.p4.fdStepDeg
%      with tau = tau*, started (continuationOptions) from the joint
%      solution so that it stays in the same local optimum.
%   3. per-spacecraft validation (validateSolution, unchanged gates) with
%      the shared formation phase, on the parallel pool.
nS = numel(p.fleet.point);
cont = continuationOptions();

% 1. Joint NLP from the best grid points
tB = tic;
fo = buildFleetOcp(p, ref, cont);
F.buildTime = toc(tB);
fprintf('Problem 4: joint NLP with %d variables and %d constraints (built in %.1f s)\n', ...
        fo.nW, fo.nG, F.buildTime);
starts = struct('phiGridDeg', {}, 'J_grid_Ns', {}, 'phiTauDeg', {}, 'tau_s', {}, 'J_Ns', {}, ...
                'status', {}, 'iterations', {}, 'time_s', {});
best = [];
for k = 1:numel(S.gridBest)
    m = S.gridBest(k);
    [sols, in] = solveFleetOcp(fo, fo.solver, fleetStart(fo, S.sub(m, :), S.phis(m)), [], p);
    starts(k) = struct('phiGridDeg', rad2deg(S.phis(m)), 'J_grid_Ns', S.J4(m), ...
        'phiTauDeg', rad2deg(in.phiTau), 'tau_s', in.tau, 'J_Ns', in.J, 'status', in.status, ...
        'iterations', in.iter, 'time_s', in.wallTime);
    fprintf('Problem 4 joint start %d/%d (grid %.1f deg, J4 %.4f N s): phiTau %.3f deg, J %.4f N s, %s (%d it., %.1f s)\n', ...
            k, numel(S.gridBest), rad2deg(S.phis(m)), S.J4(m), rad2deg(in.phiTau), in.J, ...
            in.status, in.iter, in.wallTime);
    if in.success && (isempty(best) || in.J < best.info.J)
        best = struct('sols', {sols}, 'info', in, 'k', k);
    end
end
views = fo.view;
clear fo                                % free the joint NLP before the pool
if isempty(best)
    error('LIFE:noJoint', 'No Problem 4 joint start reached Solve_Succeeded.');
end
F.starts = starts;
F.bestStart = best.k;
info = best.info; sols = best.sols;
F.jointTime = sum([starts.time_s]);

% Single-spacecraft Problem 3 OCP per pool worker (finite differences)
pool = fleetPool(p);
pc = p; pc.ocp.ipopt.print_level = 0;
for f = fieldnames(cont)', pc.ocp.ipopt.(f{1}) = cont.(f{1}); end
Cc = parallel.pool.Constant(@() buildOcp(pc, ref, 3, true));

% 2. Finite-difference check of the sensitivities: every subproblem at
% phi* + m*h, m = -2..2 (tau = tau*), started from its joint block, and the
% slope at phi* of a least-squares quadratic through the five costs.
% The re-solves form a smooth curve (rms about the fit <= 0.004 N s) but
% land up to 0.05 N s off the joint cost, also at phi* itself (nearby
% optima of the redundant thruster allocation), so central differences
% with small steps are dominated by that offset (measured: collector 1
% at +-0.1, 0.03, 0.01 deg gave -12.3, -16.1, -12.4 N s/rad, multiplier
% -16.75; at +-1 deg it changed local optimum). A point that does not
% reach Solve_Succeeded is solved again from the re-solve at phi*.
h = deg2rad(p.ocp.p4.fdStepDeg); off = h*(-2:2); nO = numel(off);
wJ = cell(1, nS); k = 0; who = zeros(nS*nO, 2);
for i = 1:nS
    v = views(i);
    blk = v.ix.X(1):v.ix.U(end);
    wJ{i} = [info.x(1); info.x(blk); 0];
    for m = 1:nO
        k = k + 1; who(k, :) = [i m];
        ws = wJ{i}; ws(end) = (info.phiTau + off(m))/pi;
        Fd(k) = parfeval(pool, @fleetWorker, 1, Cc, pc, ref, i, info.phiTau + off(m), ...
                         info.tau, ws); %#ok<AGROW>
    end
end
fd = fetchOutputs(Fd);
mc = find(off == 0);
ctr = (who(:, 1)' - 1)*nO + mc;          % re-solve at phi* of the same spacecraft
redo = find(~[fd.success] & [fd(ctr).success]);
for j = 1:numel(redo)
    k = redo(j); i = who(k, 1);
    ws = fd(ctr(k)).w; ws(end) = (info.phiTau + off(who(k, 2)))/pi;
    Fr(j) = parfeval(pool, @fleetWorker, 1, Cc, pc, ref, i, info.phiTau + off(who(k, 2)), ...
                     info.tau, ws); %#ok<AGROW>
end
if ~isempty(redo), fd(redo) = fetchOutputs(Fr); end
clear Cc
Jfd = reshape([fd.J], nO, nS)';         % nS x nO
A = [ones(nO, 1), off', off'.^2];
cf = A \ Jfd';
F.fd = struct('stepDeg', p.ocp.p4.fdStepDeg, 'offsetsDeg', rad2deg(off), 'J', Jfd, ...
    'dJdphi', cf(2, :), 'curvature', 2*cf(3, :), 'rms', sqrt(mean((A*cf - Jfd').^2, 1)), ...
    'resolveOffset', Jfd(:, mc)' - info.Ji, 'retried', numel(redo), ...
    'success', [fd.success], 'status', {cellstr([fd.status])});
fprintf(['dJ_i/dphi [N s/rad]: multipliers %s; finite differences (quadratic fit over +-%g deg) %s, ' ...
         'rms %s N s, re-solve at phi* - joint %s N s\n'], mat2str(info.dJdphi, 4), 2*p.ocp.p4.fdStepDeg, ...
        mat2str(F.fd.dJdphi, 4), mat2str(F.fd.rms, 2), mat2str(F.fd.resolveOffset, 2));

% 3. Validation of every spacecraft with the shared formation phase
[res, F.validationTime] = validateFleet(sols, info, ref, p, 4);
F.res = res;
F.info = rmfield(info, {'x', 'lamX', 'lamG'});
F.x = info.x;
end

function w0 = fleetStart(fo, subs, phi)
%FLEETSTART Joint initial point from the subproblem solutions subs{i} at
%   the phase phi (tau = tf): the blocks side by side.
w0 = zeros(fo.nW, 1);
w0(1) = 1; w0(2) = phi/pi;
for i = 1:fo.nS
    v = fo.view(i);
    w0(v.ix.X(1):v.ix.U(end)) = subs{i}.w(2:end-1);
end
end
