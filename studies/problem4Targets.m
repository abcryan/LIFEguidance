%% Problem 4 for other target directions
%  The full Problem 4 pipeline (the sweep of fleetSweep.m, then the joint
%  solve of main.m) with only the target sky direction p.target.sky
%  changed; everything else (initial sky, formations, N, d, Tmax, tf) as
%  in lifeParameters(). Outputs in results/problem4_<name>/ (the default
%  target's results/problem4/ is untouched), plus targets.json there with
%  the geometric reference phases. Finally the overview figure
%  results/problem4_target_overview.png: the initial formation and the
%  target formations at arrival (default target and every saved target
%  run) in one 3D view, as the left panel of video_fleet.mp4.
%
%  The phase phiTau is measured in each target's own pointing frame O*
%  (x_* = s x z_c), so phiTau = 0 is a convention that differs per target.
%  The reference to compare phiTau* with is the phase at which the target
%  rectangle best overlays the (co-scaled) initial one: the minimum of the
%  summed collector displacement ("phiRefDisp"), and, for information,
%  the phase of the smallest formation rotation ("phiRefRot"). For the
%  default target (both skies in the equatorial plane) both are 0.
%
%  Selected targets: set runTargets before calling (one name or several),
%  e.g.  runTargets = "A"; run studies/problem4Targets   (~45 min per
%  target); runTargets = ""; redraws only the overview figure.
targets = struct( ...
    'name',  {"A", "B"}, ...
    'thetaDeg', {80, 55}, ...           % angle from anti-Sun (+x_c at t = 0)
    'psiDeg',   {90, 225}, ...          % about x_c, from +y_c towards +z_c
    'why', {"88 deg slew towards the celestial north, direct path in the band", ...
            "48 deg slew to the south, near the inner edge (55 deg), direct path in the band"});
if ~exist('runTargets', 'var'), runTargets = [targets.name]; end

projectRoot = fileparts(fileparts(mfilename('fullpath')));   % studies/ -> project root
addpath(projectRoot);
addpath(genpath(fullfile(projectRoot, 'src')));
setupCasadi();
resultsDir = fullfile(projectRoot, 'results');

for tg = targets(ismember([targets.name], runTargets))
    tAll = tic;
    name = "problem4_target" + tg.name;
    outDir = fullfile(resultsDir, name);
    if ~isfolder(outDir), mkdir(outDir); end
    p = lifeParameters();
    p.target.sky = [cosd(tg.thetaDeg); sind(tg.thetaDeg)*cosd(tg.psiDeg); sind(tg.thetaDeg)*sind(tg.psiDeg)];
    p = deriveFormations(p);
    fprintf('\n===== Problem 4, target %s: sky [%.4f %.4f %.4f], %s =====\n', ...
            tg.name, p.target.sky, tg.why);

    % Geometry of the retargeting and the reference phases
    g = referencePhases(p);
    g.target = tg; g.sky = p.target.sky'; g.initSky = p.init.sky';
    fprintf(['Slew %.1f deg; reference phases: min collector displacement %.1f deg, ' ...
             'min formation rotation %.1f deg\n'], g.slewDeg, g.phiRefDispDeg, g.phiRefRotDeg);

    ref = propagateReference(p, p.tf);
    checkSunFeasibility(p, ref);

    % 1. Sweep and per-spacecraft Problem 3 (fleetSweep.m)
    S = runFleetSweep(p, ref);
    save(fullfile(outDir, 'sweep.mat'), 'S', '-v7.3');
    printFleetSweep(S, p);
    writeSweepJson(S, p, fullfile(outDir, 'sweep.json'));

    % 2. Joint solve, checks, validation, outputs (main.m, Problem 4)
    [F, M] = runProblem4(p, ref, resultsDir, char(name));

    % 3. Phase relative to the geometric reference
    g.phiTauDeg = rad2deg(angle(exp(1i*F.info.phiTau)));
    g.phiTauMinusRefDispDeg = rad2deg(angle(exp(1i*(F.info.phiTau - deg2rad(g.phiRefDispDeg)))));
    g.phiTauMinusRefRotDeg  = rad2deg(angle(exp(1i*(F.info.phiTau - deg2rad(g.phiRefRotDeg)))));
    g.runtime_min = toc(tAll)/60;
    fid = fopen(fullfile(outDir, 'targets.json'), 'w');
    fprintf(fid, '%s', jsonencode(g, 'PrettyPrint', true));
    fclose(fid);
    fprintf(['Target %s: phiTau* = %.2f deg, minus min-displacement reference %.2f deg, ' ...
             'minus min-rotation reference %.2f deg (%.1f min)\n'], tg.name, g.phiTauDeg, ...
            g.phiTauMinusRefDispDeg, g.phiTauMinusRefRotDeg, g.runtime_min);
    close all
end

% Close the pool explicitly (MATLAB crashed shutting it down at exit)
delete(gcp('nocreate'));

% Overview of all saved runs
plotTargetOverview(resultsDir, [targets.name]);

function g = referencePhases(p)
% Phases (in the target's O*) at which the target formation at arrival is
% closest to the initial formation at t = 0: the summed collector
% displacement (the initial pyramid is the target one scaled, so the
% scale does not move the minimum) and the rotation angle between them.
Ri = skyBasis(p.init.sky)*rot3(p.init.phase0);
Rt = skyBasis(p.target.sky);
Pt = pyramidPoints(p.target); Pi = pyramidPoints(p.init);
ph = deg2rad(-180:0.05:180);
dsp = zeros(size(ph)); rot = dsp;
for k = 1:numel(ph)
    Rk = Rt*rot3(ph(k));
    dsp(k) = sum(vecnorm(Rk*Pt(:, 1:4) - Ri*Pi(:, 1:4)));
    rot(k) = acos(min(1, (trace(Ri'*Rk) - 1)/2));
end
[~, i] = min(dsp); [~, j] = min(rot);
g.slewDeg = acosd(p.init.sky'*p.target.sky);
g.phiRefDispDeg = rad2deg(ph(i));
g.phiRefRotDeg = rad2deg(ph(j));
end

function plotTargetOverview(resultsDir, names)
% Initial formation at t = 0 (grey) and the target formations at arrival
% of the default target (results/problem4) and of every saved target run,
% SER-parallel axes centred on the formation; with each target's slew,
% optimal phase and fleet impulse.
dirs = ["problem4", "problem4_target" + names];
labels = ["original", names];
col = [0.80 0.10 0.10; 0.05 0.45 0.80; 0.10 0.60 0.30; 0.85 0.45 0.05; 0.45 0.20 0.60];
fig = figure('Color', 'w', 'Position', [80 80 1250 950]);
ax = axes(fig); hold(ax, 'on');
hl = gobjects(0); P = [];
for k = 1:numel(dirs)
    f = fullfile(resultsDir, dirs(k), 'solution.mat');
    if ~isfile(f), fprintf('Overview: %s missing, skipped.\n', f); continue, end
    d = load(f, 'fleet');
    pk = d.fleet.res{1}.p;              % its p carries the solution's phase
    tau = d.fleet.info.tau;
    if isempty(hl)                      % initial formation, common to all runs
        V0 = initialFormationPoints(pk);
        h0 = drawFormation(ax, V0, [0.55 0.55 0.55]);
        R0 = frames(0, pk).Rbf;
        setFormationPose(h0, R0);
        X0 = R0*V0; P = [P, X0]; %#ok<AGROW>
        h0.dots.DisplayName = 'initial formation, t = 0';
        hl(end + 1) = h0.dots; %#ok<AGROW>
        labelPoints(ax, X0, [0.45 0.45 0.45], pk.fleet.point);
    end
    V = [pk.viz.collectors, pk.viz.combiner];
    h = drawFormation(ax, V, col(k,:));
    Rf = frames(tau, pk).Rbf;
    setFormationPose(h, Rf);
    X = Rf*V; P = [P, X]; %#ok<AGROW>
    h.dots.DisplayName = sprintf('target %s at t = %.0f h: slew %.0f%s, \\phi_\\tau^* = %.1f%s, J_4^* = %.0f N s', ...
        labels(k), tau/3600, acosd(pk.init.sky'*pk.target.sky), char(176), ...
        rad2deg(angle(exp(1i*d.fleet.info.phiTau))), char(176), d.fleet.info.J);
    hl(end + 1) = h.dots; %#ok<AGROW>
    labelPoints(ax, X, col(k,:), pk.fleet.point);
end
if isempty(hl), close(fig); return, end
plot3(ax, 0, 0, 0, 'k+', 'MarkerSize', 12, 'LineWidth', 1.5);
axis(ax, 'equal'); axis(ax, zoomLimits(P, 60, 0.15)); grid(ax, 'on'); box(ax, 'on');
view(ax, 35, 20);
xlabel(ax, 'x_b (m)'); ylabel(ax, 'y_b (m)'); zlabel(ax, 'z_b (m)');
title(ax, {'Problem 4: initial formation and the target formations at arrival', ...
    'SER-parallel axes, origin at the formation centre (+); numbers: collectors 1-4, combiner 5 (apex)'}, ...
    'FontWeight', 'normal');
set(ax, 'FontSize', 11);
legend(ax, hl, 'Location', 'southoutside', 'FontSize', 11, 'Box', 'off');
f = fullfile(resultsDir, 'problem4_target_overview.png');
exportgraphics(fig, f, 'Resolution', 200);
fprintf('Written %s\n', f);
end

function labelPoints(ax, X, c, pt)
% Spacecraft numbers (document numbering; pt = p.fleet.point maps
% spacecraft i to its pyramid point) next to the points.
for i = 1:numel(pt)
    q = X(:, pt(i));
    text(ax, q(1), q(2), q(3), sprintf('  %d', i), 'Color', c, 'FontWeight', 'bold', 'FontSize', 10);
end
end
