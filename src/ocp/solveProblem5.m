function F = solveProblem5(p, ref, P4)
%SOLVEPROBLEM5 Problem 5 (Secs. 4.6, 5.2): Problem 4 with collision
%   avoidance |r_i - r_j| >= dmin for all ten pairs at every collocation
%   point (imposed with dmin + margin), from the Problem 4 optimum P4
%   (results/problem4/solution.mat: fleet, metrics). One compiled solver
%   (buildFleetOcp with collision rows) serves every dmin, which only sets
%   the bounds (setCollision).
%   1. Starts at dmin (tau and phiTau free, primal start with
%      continuationOptions as in Problem 4): the Problem 4 optimum, and
%      the starts in which the two collectors of each crossing pair
%      (Problem 4 separation below dmin + margin) are displaced in opposite
%      directions along +-z_* by (dmin/2) sin(pi nu), for every
%      combination of passing sides. Best Solve_Succeeded kept. Fallback
%      if none converges: continuation dmin = 5, 10, 15, ... from Problem 4.
%   2. The sweep J5(dmin) on the chosen branch, continued from the dmin
%      solution downwards and upwards (p.collision.sweep), and the central
%      differences at dmin +- fdStep for dJ5/ddmin.
%   3. dJ5/ddmin from the collision multipliers (envelope theorem).
%   4. Validation of every spacecraft (validateFleet, unchanged gates).
%   Serial in the client (p.collision.nWorkers = 0, one NLP), or with the
%   starts and the sweep chains in parallel on a pool of
%   p.collision.nWorkers workers, each with its own compiled NLP (~4 GB
%   each at N = 400, d = 5), with the same solves and results.
cp = p.collision; dmin = cp.dmin; mg = cp.margin;
cont = continuationOptions();
nWk = 0; if isfield(cp, 'nWorkers'), nWk = cp.nWorkers; end

% Crossing pairs: closer than dmin + margin in the Problem 4 solution
d4 = P4.metrics.dist;
[d4min, k4] = min(d4.D, [], 2);
cross = find(d4min < dmin + mg)';
F.crossing = struct('pair', num2cell(d4.pairs(cross, :), 2)', 'p4MinM', num2cell(d4min(cross))', ...
                    'p4MinTime', num2cell(d4.t(k4(cross))));
fprintf('Crossing pairs (Problem 4 separation below %.1f m): %s\n', dmin + mg, ...
        strjoin(arrayfun(@(k) sprintf('(%d,%d) %.2f m', d4.pairs(k,:), d4min(k)), cross, ...
        'UniformOutput', false), ', '));

% CHOICE (p.collision.pairs): "all" imposes (67) on all ten pairs as in the
% document; "crossing" only on the crossing pairs (4000 instead of 20000
% rows), the other pairs (>= 128 m apart in Problem 4) are verified on the
% propagated trajectories only (problem5Metrics, step (iii), all pairs).
colPairs = true;
if isfield(cp, 'pairs') && strcmp(cp.pairs, 'crossing') && ~isempty(cross)
    colPairs = d4.pairs(cross, :);
end

% Problem 5 NLP: in the client, or once per pool worker
tB = tic;
if nWk == 0
    fo = buildFleetOcp(p, ref, cont, colPairs);
    meta = struct('pairs', fo.pairs, 'iCol', fo.iCol, 'colPair', fo.colPair, ...
                  'colTime', fo.colTime, 'nW', fo.nW, 'nG', fo.nG);
    run = @(jobs) cellfun(@(a) problem5Solve(fo, p, a{:}), jobs, 'UniformOutput', false);
else
    pool = gcp('nocreate');
    if ~isempty(pool) && pool.NumWorkers ~= nWk, delete(pool); end
    pc = p; pc.par.nWorkers = nWk; pc.ocp.ipopt.print_level = 0;
    pool = fleetPool(pc);
    C = parallel.pool.Constant(@() buildFleetOcp(pc, ref, cont, colPairs));
    m = fetchOutputs(parfevalOnAll(pool, @problem5Worker, 1, C, 'meta'));
    meta = m(1);
    run = @(jobs) runParallel(pool, C, p, jobs);
end
F.buildTime = toc(tB);
fprintf(['Problem 5: joint NLP with %d variables, %d constraints, %d collision rows for %d pairs ' ...
         '(built in %.1f s, %s)\n'], meta.nW, meta.nG, numel(meta.iCol), size(meta.pairs, 1), ...
        F.buildTime, modeLabel(nWk));
F.pairs = meta.pairs;

% 1. Starts at dmin
w4 = P4.fleet.x;
nc = numel(cross);                       % every combination of passing sides
sides = [zeros(1, nc); 1 - 2*(dec2bin(0:2^nc - 1, max(nc, 1)) - '0')];
sides = sides(:, 1:nc);
nStart = size(sides, 1);
names = cell(1, nStart); jobs = cell(1, nStart);
for k = 1:nStart
    if k == 1
        names{k} = 'Problem 4 optimum'; jobs{k} = {w4, dmin};
    else
        names{k} = ['displaced ' strjoin(arrayfun(@(m) sprintf('(%d,%d) %+d', d4.pairs(cross(m), :), ...
                    sides(k, m)), 1:nc, 'UniformOutput', false), ', ')];
        jobs{k} = {w4, dmin, sides(k, :), d4.pairs(cross, :)};
    end
end
R = run(jobs);
starts = struct('name', {}, 'sides', {}, 'J_Ns', {}, 'tau_s', {}, 'phiTauDeg', {}, 'status', {}, ...
                'iterations', {}, 'time_s', {}, 'minSepM', {});
best = [];
for k = 1:nStart
    r = R{k}; in = r.info;
    starts(k) = struct('name', names{k}, 'sides', sides(k, :), 'J_Ns', in.J, 'tau_s', in.tau, ...
        'phiTauDeg', rad2deg(in.phiTau), 'status', in.status, 'iterations', in.iter, ...
        'time_s', in.wallTime, 'minSepM', min(r.cs.minPerPair));
    fprintf('Problem 5 start %d/%d (%s): J %.4f N s, phiTau %.3f deg, min sep %.3f m, %s (%d it., %.1f s)\n', ...
            k, nStart, names{k}, in.J, rad2deg(in.phiTau), starts(k).minSepM, in.status, in.iter, in.wallTime);
    if in.success && (isempty(best) || in.J < best.info.J)
        best = r; best.k = k;
    end
end
F.continuationFallback = [];
if isempty(best)
    % Fallback: continuation in dmin from the Problem 4 optimum
    steps = unique([5:5:dmin, dmin]);
    fprintf('Problem 5: no start converged; continuation in dmin %s\n', mat2str(steps));
    Rc = run({{w4, steps}}); Rc = Rc{1};
    for k = 1:numel(Rc)
        F.continuationFallback(k).dmin = Rc(k).dmin;
        F.continuationFallback(k).status = Rc(k).info.status;
        F.continuationFallback(k).J_Ns = Rc(k).info.J;
        fprintf('  dmin %4.1f m: J %.4f N s, %s (%d it., %.1f s)\n', Rc(k).dmin, Rc(k).info.J, ...
                Rc(k).info.status, Rc(k).info.iter, Rc(k).info.wallTime);
    end
    if ~Rc(end).info.success
        error('LIFE:noProblem5', 'Problem 5 did not converge from any start or the continuation.');
    end
    best = Rc(end); best.k = 0;
end
F.starts = starts; F.bestStart = best.k;
info = best.info; sols = best.sols;
fprintf('Problem 5: best start %d, J5 = %.4f N s\n', best.k, info.J);

% 2. Sweep J5(dmin) on this branch (two chains from the dmin solution) and
% central differences at dmin +- fdStep
sw = sort(cp.sweep);
down = fliplr(sw(sw < dmin)); up = sw(sw > dmin);
fdD = dmin + cp.fdStep*[-1 1];
% plus the smallest sweep distance from the Problem 4 optimum: there the
% constraint is inactive, so the NLP must return J4* (direct check of the
% equivalence, independent of the branch the downward chain follows)
dLow = sw(1);
jobs = {{info.x, down}, {info.x, up}, {info.x, fdD(1)}, {info.x, fdD(2)}, {w4, dLow}};
use = [~isempty(down), ~isempty(up), true, true, true];
R2 = cell(1, 5); R2(use) = run(jobs(use));
F.sweep = sweepRow(best);
for c = find(use(1:2))
    for r = R2{c}
        F.sweep(end + 1) = sweepRow(r);
        fprintf('Problem 5 sweep dmin %4.1f m: J %.4f N s, dJ/ddmin %.4f N s/m, %d active pairs, %s (%d it., %.1f s)\n', ...
                r.dmin, r.info.J, r.sens.dJddmin, numel(r.sens.activePairs), r.info.status, ...
                r.info.iter, r.info.wallTime);
    end
end
[~, o] = sort([F.sweep.dmin]); F.sweep = F.sweep(o);
Jfd = [R2{3}.info.J, R2{4}.info.J]; okFd = [R2{3}.info.success, R2{4}.info.success];
for m = 1:2
    in = R2{2 + m}.info;
    fprintf('Problem 5 finite difference dmin %.2f m: J %.4f N s, %s (%d it., %.1f s)\n', ...
            fdD(m), in.J, in.status, in.iter, in.wallTime);
end
F.fd = struct('dmin', fdD, 'J', Jfd, 'success', okFd, 'dJddmin', diff(Jfd)/diff(fdD));
in = R2{5}.info;
F.lowFromP4 = struct('dmin', dLow, 'J_Ns', in.J, 'success', in.success, 'status', in.status, ...
                     'iterations', in.iter, 'nActivePairs', numel(R2{5}.sens.activePairs));
fprintf('Problem 5 at dmin %.1f m from the Problem 4 optimum: J %.4f N s (J4* %.4f), %s (%d it., %.1f s)\n', ...
        dLow, in.J, P4.fleet.info.J, in.status, in.iter, in.wallTime);

% 3. Sensitivity and active set at the optimum
F.sens = best.sens;
fprintf('dJ5/ddmin: multipliers %.4f N s/m, central differences %.4f N s/m\n', ...
        F.sens.dJddmin, F.fd.dJddmin);
F.colloc = best.cs;

% 4. Validation of every spacecraft (the workers' NLPs are released first)
if nWk > 0, clear C; end
[F.res, F.validationTime] = validateFleet(sols, info, ref, p, 5);
F.info = rmfield(info, {'x', 'lamX', 'lamG'});
F.x = info.x;
F.lamCol = info.lamG(meta.iCol);
F.colPair = meta.colPair; F.colTime = meta.colTime;
F.nWorkers = nWk;
F.problem = 5;
end

function R = runParallel(pool, C, p, jobs)
% Every job on the pool (problem5Worker), results in job order.
for k = numel(jobs):-1:1
    fut(k) = parfeval(pool, @problem5Worker, 1, C, 'solve', p, jobs{k}{:});
end
R = cell(1, numel(jobs));
for k = 1:numel(jobs)
    R{k} = fetchOutputs(fut(k));
end
end

function r = sweepRow(s)
r = struct('dmin', s.dmin, 'J_Ns', s.info.J, 'success', s.info.success, 'status', s.info.status, ...
    'iterations', s.info.iter, 'time_s', s.info.wallTime, 'tau_s', s.info.tau, ...
    'phiTauDeg', rad2deg(s.info.phiTau), 'dJddmin', s.sens.dJddmin, ...
    'nActivePairs', numel(s.sens.activePairs), 'minSepM', min(s.cs.minPerPair));
end

function s = modeLabel(n)
if n == 0, s = 'serial'; else, s = sprintf('on each of %d pool workers', n); end
end
