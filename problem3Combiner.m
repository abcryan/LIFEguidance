%% Problem 3 for the combiner alone (spacecraft 5), vs collector 1
%  The combiner (apex of the pyramid, 1246 m from its target point) solved
%  as Problem 3 with the same procedure as collector 1 in main.m: the
%  8-phase multi-start from the common initial guess (Sec. 5.2), the
%  best Solve_Succeeded start kept. If results/problem4/sweep.mat matches
%  the settings, one more start from the combiner's best sweep point
%  (tau and phi freed, as in fleetSweep.m) is added and reported
%  separately. Validation, outputs and figures as for Problems 1-3, in
%  results/problem3_combiner/, plus a comparison with collector 1's
%  Problem 3 (results/problem3, run main.m first).
projectRoot = fileparts(mfilename('fullpath'));
addpath(projectRoot);
addpath(genpath(fullfile(projectRoot, 'src')));
setupCasadi();

iSc = 5;                                % the combiner
p = lifeParameters();
ref = propagateReference(p, p.tf);
pk = spacecraftParams(p, iSc);
checkSunFeasibility(pk, ref);
resultsDir = fullfile(projectRoot, 'results');
outDir = fullfile(resultsDir, 'problem3_combiner');

% 1. The multi-start of main.m
ocp = buildOcp(pk, ref, 3);
[sol, info, starts] = solveMultiStart(ocp, pk);
info.source = 'multi-start';
msBest = struct('J_Ns', sol.J, 'phiTauDeg', rad2deg(sol.phiTau));

% 2. One start from the combiner's best sweep point (fleetSweep.m)
sweepFile = fullfile(resultsDir, 'problem4', 'sweep.mat');
S = [];
if isfile(sweepFile)
    d = load(sweepFile, 'S');
    if strcmp(d.S.key, fleetKey(p)), S = d.S; end
end
warmStart = [];
if ~isempty(S)
    Ji = S.Jbest(:, iSc); Ji(S.isPresc) = NaN;
    [~, mb] = min(Ji);
    sb = S.sub{mb, iSc};
    pw = pk;
    wo = warmStartOptions();
    for f = fieldnames(wo)', pw.ocp.ipopt.(f{1}) = wo.(f{1}); end
    clear ocp
    ocpW = buildOcp(pw, ref, 3);
    lam = struct('lam_x0', [sb.lamX(1:end-1); 0], 'lam_g0', sb.lamG);   % phase now free
    [sW, inW] = solveOcp(ocpW, sb.w, pw, [], lam);   % same w layout as the parametric OCP
    warmStart = struct('fromDeg', rad2deg(S.phis(mb)), 'phiTauDeg', rad2deg(sW.phiTau), ...
        'J_Ns', sW.J, 'tau_s', sW.tau, 'status', inW.status, 'iterations', inW.iter, ...
        'time_s', inW.wallTime);
    fprintf('Problem 3 start from best sweep point %.1f deg: phiTau -> %.2f deg, J %.4f N s, %s (%d it., %.1f s)\n', ...
            warmStart.fromDeg, warmStart.phiTauDeg, sW.J, inW.status, inW.iter, inW.wallTime);
    if inW.success && sW.J < sol.J
        sol = sW; inW.multiStartTime = info.multiStartTime + inW.wallTime;
        inW.bestStart = 0; inW.source = 'best sweep point'; info = inW;
    end
else
    fprintf('No matching results/problem4/sweep.mat: multi-start only.\n');
end

% 3. Validation and outputs, as main.m
pv = pk;
pv.target.phase0 = sol.phase0; pv.phi0 = sol.phase0;
res = validateSolution(sol, info, ref, pv, 3);
res.multistart = starts;
res.warmStart = warmStart;
res.multistartBest = msBest;
res.cls = pathClass(res.sun.azimuth);
saveResults(res, ref, pv, outDir);
plotResults(res, pv, outDir);
printSummary(res, pv);

% 4. Comparison with collector 1 (Problem 3 of main.m)
f1 = fullfile(resultsDir, 'problem3', 'solution.mat');
if ~isfile(f1)
    fprintf('No %s: run main.m for the comparison.\n', f1);
    return
end
c1 = load(f1, 'res', 'p');
R = {c1.res, res}; P = {c1.p, pv};
names = ["collector 1 (+,+)", "combiner"];
for k = 2:-1:1
    r = R{k}; q = P{k}; s = r.sol;
    x0 = initialRelativeState(q);
    qa = s.X(1:4, 1); qb = s.X(1:4, end);
    row(k) = struct('spacecraft', names(k), ...
        'initial_distance_m', norm(r.prop.x(8:10, 1)), ...
        'J_Ns', r.J, 'propellant_kg', r.dm, ...
        'multistart_best_J_Ns', NaN, ...
        'phi_tau_deg', rad2deg(angle(exp(1i*s.phiTau))), 'tau_s', s.tau, ...
        'tau_at_upper_bound', r.info.tauAtUpper, ...
        'attitude_change_deg', rad2deg(2*acos(min(1, abs(qa'*qb)))), ...
        'max_rate_deg_s', rad2deg(max(vecnorm(s.X(5:7, :)))), ...
        'max_thrust_N', r.maxThrust, 'n_burns', numel(r.burns), ...
        'path_class', pathClass(r.sun.azimuth), ...
        'anti_sun_min_deg', min(r.sun.antiDeg), 'anti_sun_max_deg', max(r.sun.antiDeg), ...
        'sun_active_points', numel(r.sun.active), ...
        'iterations', r.info.iter, 'multistart_time_s', r.info.multiStartTime, ...
        'accepted', r.accepted);
end
row(1).multistart_best_J_Ns = c1.res.J;
row(2).multistart_best_J_Ns = msBest.J_Ns;
cmp.table = struct2table(row);
cmp.ratio_J = row(2).J_Ns/row(1).J_Ns;
cmp.ratio_distance = row(2).initial_distance_m/row(1).initial_distance_m;
if ~isempty(S)
    cmp.sweep = struct('phiDeg', rad2deg(S.phis(~S.isPresc)), ...
        'J_collector1_Ns', S.Jbest(~S.isPresc, 1)', 'J_combiner_Ns', S.Jbest(~S.isPresc, iSc)');
    cmp.collector1_best_P3_fleetSweep_Ns = S.p3(1).J;
end
fid = fopen(fullfile(outDir, 'comparison_collector1.json'), 'w');
fprintf(fid, '%s', jsonencode(cmp, 'PrettyPrint', true));
fclose(fid);
fprintf('\nProblem 3, collector 1 vs combiner\n');
disp(rows2vars(cmp.table(:, 2:end), 'VariableNamingRule', 'preserve'));
fprintf('Combiner / collector 1: J x%.2f, initial distance x%.2f\n', cmp.ratio_J, cmp.ratio_distance);

fig = plotCombinerComparison(R, P, S, iSc);
exportgraphics(fig, fullfile(outDir, 'comparison_collector1.png'), 'Resolution', 200);

function fig = plotCombinerComparison(R, P, S, iSc)
% 2x2: cumulative impulse, J(phiTau) - min from the fleet sweep, anti-Sun
% angle with the band, distance to the target point (log).
st = plotStyle();
col = [st.sc(1,:); st.sc(iSc,:)];
lab = ["collector 1", "combiner"];
fig = figure('Color', 'w', 'Position', [100 100 1300 850]);
tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, 'Problem 3: collector 1 vs combiner (spacecraft 5)');

ax = nexttile(tl); hold(ax, 'on'); grid(ax, 'on');
for k = 1:2
    s = R{k}.sol;
    tb = [0, cumsum(s.h)];
    plot(ax, tb/3600, R{k}.cumImpulse, 'Color', col(k,:), 'LineWidth', 1.6, ...
         'DisplayName', sprintf('%s: %.1f N s', lab(k), R{k}.J));
end
xlabel(ax, 't [h]'); ylabel(ax, 'cumulative impulse [N s]');
title(ax, '(a) impulse'); legend(ax, 'Location', 'northwest');

ax = nexttile(tl); hold(ax, 'on'); grid(ax, 'on');
if ~isempty(S)
    ph = rad2deg(S.phis(~S.isPresc));
    for k = 1:2
        i = [1 iSc]; J = S.Jbest(~S.isPresc, i(k))';
        plot(ax, ph, J - min(J), '-o', 'Color', col(k,:), 'LineWidth', 1.4, ...
             'MarkerFaceColor', col(k,:), 'MarkerSize', 4, ...
             'DisplayName', sprintf('%s (min %.1f N s)', lab(k), min(J)));
    end
    xlim(ax, [-180 180]); xticks(ax, -180:45:180);
    legend(ax, 'Location', 'north');
else
    text(ax, 0.5, 0.5, 'no fleet sweep', 'HorizontalAlignment', 'center');
end
xlabel(ax, '\phi_\tau [deg]'); ylabel(ax, 'J(t_f, \phi_\tau) - min J [N s]');
title(ax, '(b) phase dependence (fixed-phase sweep of fleetSweep.m)');

ax = nexttile(tl); hold(ax, 'on'); grid(ax, 'on');
q = P{1}; band = sort(180 - rad2deg([q.sun.thetaMin, q.sun.thetaMax]));
patch(ax, [0 4 4 0], band([1 1 2 2]), st.band, 'FaceAlpha', 0.15, 'EdgeColor', 'none', ...
      'DisplayName', 'band');
for k = 1:2
    plot(ax, R{k}.sun.t/3600, R{k}.sun.antiDeg, 'Color', col(k,:), 'LineWidth', 1.4, ...
         'DisplayName', lab(k));
end
xlabel(ax, 't [h]'); ylabel(ax, 'anti-Sun angle [deg]');
title(ax, '(c) Sun angle'); legend(ax, 'Location', 'best');

ax = nexttile(tl); hold(ax, 'on'); grid(ax, 'on');
for k = 1:2
    pr = R{k}.prop;
    semilogy(ax, pr.t/3600, max(vecnorm(pr.x(8:10,:)), 1e-6), 'Color', col(k,:), ...
             'LineWidth', 1.4, 'DisplayName', lab(k));
end
set(ax, 'YScale', 'log');
xlabel(ax, 't [h]'); ylabel(ax, '|r^f_{pi}| [m]');
title(ax, '(d) distance to the target point'); legend(ax, 'Location', 'southwest');
end
