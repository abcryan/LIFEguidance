function S = runFleetSweep(p, ref)
%RUNFLEETSWEEP Problem 4 initialization and bounds (Sec. 5.2), on a
%   parallel pool (fleetPool) with one compiled parametric Problem 3 OCP
%   per worker:
%   1. Sweep: the five single-spacecraft subproblems at tau = tf and each
%      fixed phase of p.ocp.p4.sweepPhases, plus Problem 2's prescribed
%      phase phi0 + rate*tf (an upper bound for Problem 4), each from the
%      two attitude-class guesses (initialGuess, via = +-1); the better
%      one is kept. Where both fail, two retries: the direct-slerp guess
%      and a warm start from the nearest successful phase of the same
%      spacecraft. J4(phi) = sum_i J_i(tf, phi), eq. (65).
%   2. Problem 3 for every spacecraft (J*_3,i): the 8-phase multi-start
%      of Problem 3 from the common guess, plus one start from the
%      spacecraft's best sweep point (tau and phi freed, primal-dual warm
%      start, so J*_3,i <= min_phi J_i(tf, phi)); the best is kept.
%   S holds every solve (costs, status, classes, the sweep sensitivities
%   dJ_i/dphi from the phase multiplier) and, for the kept solve of each
%   (phi, i), the solution and multipliers (joint warm start).
p.ocp.ipopt.print_level = 0;
pool = fleetPool(p);
fprintf('Problem 4 sweep on %d workers: compiling the Problem 3 OCP on each ...\n', pool.NumWorkers);
C = parallel.pool.Constant(@() buildOcp(p, ref, 3, true));
nS = numel(p.fleet.point);
S.phiPresc = p.target.phase0 + p.target.rate*p.tf;
S.phis = [p.ocp.p4.sweepPhases, S.phiPresc];
S.isPresc = [false(1, numel(p.ocp.p4.sweepPhases)), true];
nP = numel(S.phis); vias = [1 -1];
S.vias = vias;
starts = p.ocp.p3Starts; nM = numel(starts);
tStart = tic;
run = @(args, labels) collect(arrayfun(@(k) parfeval(pool, @fleetWorker, 1, C, args{k}{:}), ...
                                       1:numel(args)), labels);

% 1. Sweep, and the Problem 3 multi-start (independent) in the same queue.
% Task order: guess fastest, then spacecraft, then phase.
args = {}; labels = {};
for m = 1:nP
    for i = 1:nS
        for v = 1:2
            args{end + 1} = {p, ref, i, S.phis(m), p.tf, struct('phi', S.phis(m), 'via', vias(v))}; %#ok<AGROW>
            labels{end + 1} = sprintf('sc %d sweep phi %7.1f deg, via %+d', i, rad2deg(S.phis(m)), vias(v)); %#ok<AGROW>
        end
    end
end
nSweep = numel(args);
for i = 1:nS
    for m = 1:nM
        args{end + 1} = {p, ref, i, NaN, NaN, struct('phi', starts(m), 'via', 0)}; %#ok<AGROW>
        labels{end + 1} = sprintf('sc %d P3 start %7.1f deg', i, rad2deg(starts(m))); %#ok<AGROW>
    end
end
R = run(args, labels);
sw = R(1:nSweep); ms = R(nSweep + 1:end);

perm = @(f) permute(reshape([sw.(f)], 2, nS, nP), [3 2 1]);   % -> nP x nS x 2
S.J = perm('J'); S.success = perm('success'); S.iter = perm('iter');
S.time = perm('time'); S.cls = perm('cls'); S.dJdphi = perm('dJdphi');
S.status = perm('status');

% Kept solve per (phi, i): the lower cost among the successful guesses
Jok = S.J; Jok(~S.success) = inf;
[Jb, S.vBest] = min(Jok, [], 3);
S.sub = cell(nP, nS); S.Jbest = NaN(nP, nS);
S.clsBest = zeros(nP, nS); S.dJbest = NaN(nP, nS);
for m = 1:nP
    for i = 1:nS
        if isinf(Jb(m, i)), continue, end
        keep(sw((m - 1)*nS*2 + (i - 1)*2 + S.vBest(m, i)), m, i);
    end
end

% Retries where both guesses failed (vBest = 3: direct slerp, 4: warm)
[mf, jf] = find(isnan(S.Jbest));
S.retry = struct('sc', {}, 'phiDeg', {}, 'kind', {}, 'J', {}, 'status', {});
if ~isempty(mf)
    args = {}; labels = {}; who = zeros(0, 3);
    for k = 1:numel(mf)
        [m, i] = deal(mf(k), jf(k));
        args{end + 1} = {p, ref, i, S.phis(m), p.tf, struct('phi', S.phis(m), 'via', 0)}; %#ok<AGROW>
        labels{end + 1} = sprintf('sc %d retry phi %7.1f deg, direct guess', i, rad2deg(S.phis(m))); %#ok<AGROW>
        who(end + 1, :) = [m i 3]; %#ok<AGROW>
        ok = find(~isnan(S.Jbest(:, i)))';
        if ~isempty(ok)
            [~, b] = min(abs(angle(exp(1i*(S.phis(ok) - S.phis(m))))));
            args{end + 1} = {p, ref, i, S.phis(m), p.tf, S.sub{ok(b), i}.w}; %#ok<AGROW>
            labels{end + 1} = sprintf('sc %d retry phi %7.1f deg, warm from %.1f deg', i, ...
                rad2deg(S.phis(m)), rad2deg(S.phis(ok(b)))); %#ok<AGROW>
            who(end + 1, :) = [m i 4]; %#ok<AGROW>
        end
    end
    Q = run(args, labels);
    for k = 1:numel(Q)
        [m, i, kind] = deal(who(k, 1), who(k, 2), who(k, 3));
        S.retry(end + 1) = struct('sc', i, 'phiDeg', rad2deg(S.phis(m)), 'kind', kind, ...
                                  'J', Q(k).J, 'status', Q(k).status);
        if Q(k).success && ~(Q(k).J >= S.Jbest(m, i))
            S.vBest(m, i) = kind;
            keep(Q(k), m, i);
        end
    end
end
S.J4 = sum(S.Jbest, 2)';                % NaN if a subproblem failed throughout
J4g = S.J4; J4g(S.isPresc) = NaN;
[~, order] = sort(J4g);
order = order(~isnan(J4g(order)));
S.gridBest = order(1:min(p.ocp.p4.nJointStarts, numel(order)));

% 2. Problem 3 per spacecraft: the multi-start, then one start from the
% best sweep point of the spacecraft, with the same OCP built with
% warm-start options (the first one is released first, to save memory)
clear C run                              % run holds C as well
pw = p;
for f = fieldnames(warmStartOptions())', pw.ocp.ipopt.(f{1}) = warmStartOptions().(f{1}); end
Cw = parallel.pool.Constant(@() buildOcp(pw, ref, 3, true));
args = {}; labels = {};
for i = 1:nS
    Ji = S.Jbest(:, i); Ji(S.isPresc) = NaN;
    [~, mb] = min(Ji);
    S.p3WarmFrom(i) = S.phis(mb);
    sb = S.sub{mb, i};
    lam = struct('lam_x0', [sb.lamX(1:end-1); 0], 'lam_g0', sb.lamG);   % phase now free
    args{end + 1} = {pw, ref, i, NaN, NaN, sb.w, lam}; %#ok<AGROW>
    labels{end + 1} = sprintf('sc %d P3 from best sweep point %.1f deg', i, rad2deg(S.phis(mb))); %#ok<AGROW>
end
W = collect(arrayfun(@(k) parfeval(pool, @fleetWorker, 1, Cw, args{k}{:}), 1:numel(args)), labels);
for i = 1:nS
    c = ms((i - 1)*nM + (1:nM));
    cand = [c, W(i)];
    st = struct('start', [cellstr(compose('phi0 = %.0f deg', rad2deg(starts))), ...
                          {sprintf('best sweep point (%.1f deg)', rad2deg(S.p3WarmFrom(i)))}], ...
        'phiTauDeg', num2cell(rad2deg([cand.phiTau])), 'J_Ns', {cand.J}, ...
        'tau_s', {cand.tau}, 'cls', {cand.cls}, 'status', cellstr([cand.status]), ...
        'iterations', {cand.iter}, 'time_s', {cand.time});
    Jc = [cand.J]; Jc(~[cand.success]) = inf;
    [~, b] = min(Jc);
    if isinf(Jc(b))
        error('LIFE:noP3', 'No Problem 3 solve of spacecraft %d succeeded.', i);
    end
    Jm = Jc(1:nM);
    S.p3(i) = struct('J', cand(b).J, 'phiTau', cand(b).phiTau, 'tau', cand(b).tau, ...
        'cls', cand(b).cls, 'best', st(b).start, 'multistartBestJ', min(Jm), ...
        'starts', {st}, 'w', cand(b).w);
end
S.wallTime = toc(tStart);
S.cpuTime = sum([R.time, W.time]);
if ~isempty(S.retry), S.cpuTime = S.cpuTime + sum([Q.time]); end
S.key = fleetKey(p);
fprintf('Problem 4 sweep: %d solves, %.1f min wall time (%.1f min solver time)\n', ...
        numel(R) + numel(W) + numel(S.retry), S.wallTime/60, S.cpuTime/60);

    function keep(r, m, i)
        S.sub{m, i} = struct('w', r.w, 'lamX', r.lamX, 'lamG', r.lamG);
        S.Jbest(m, i) = r.J;
        S.clsBest(m, i) = r.cls;
        S.dJbest(m, i) = r.dJdphi;
    end
end

function R = collect(F, labels)
%COLLECT Fetch the futures F as they finish, with a progress line each;
%   R(k) is the result of F(k).
n = numel(F); t0 = tic;
for c = 1:n
    [k, r] = fetchNext(F);
    R(k) = r; %#ok<AGROW>
    fprintf('[%3d/%3d %5.0f s] %s -> phi %7.2f deg: J %9.4f N s, class %+d, %s (%d it., %.1f s)\n', ...
            c, n, toc(t0), labels{k}, rad2deg(r.phiTau), r.J, r.cls, r.status, r.iter, r.time);
end
end
