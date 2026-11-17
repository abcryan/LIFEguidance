function figs = plotFleet(F, M, S, p, outDir)
%PLOTFLEET Problem 4 figures (Sec. 5.3), one colour per spacecraft
%   (plotStyle s.sc), saved as PNG to outDir and left open:
%   fleet_phase_sweep   J_i(tf, phi) stacked, J4(phi), phi*, phi*_3,i,
%                       per-spacecraft pull and the sensitivities
%   fleet_impulse       J_4,i vs J*_3,i, Pi4, cumulative impulse
%   fleet_trajectories  paths in the target formation frame (3D) and
%                       close-ups of every spacecraft about its target
%   fleet_tracking      norms of the target-relative errors with the gates
%   fleet_sun           anti-Sun angle and sky plot of all spacecraft
%   fleet_thrust        T_ik/Tmax of every spacecraft over time
%   fleet_distance      distances between all pairs (not constrained)
list = {
    'fleet_phase_sweep',  @() sweepFig(F, M, S, p)
    'fleet_impulse',      @() impulseFig(F, M, p)
    'fleet_trajectories', @() trajectoryFig(F, p)
    'fleet_tracking',     @() trackingFig(F, p)
    'fleet_sun',          @() sunFig(F, M, p)
    'fleet_thrust',       @() thrustFig(F, p)
    'fleet_distance',     @() distanceFig(M, p)
    };
for k = 1:size(list, 1)
    figs.(list{k, 1}) = list{k, 2}();
    exportgraphics(figs.(list{k, 1}), fullfile(outDir, [list{k, 1} '.png']), 'Resolution', 200);
end
end

function lbl = scLabel(p, i)
lbl = sprintf('%d %s', i, p.fleet.name(i));
end

function fig = sweepFig(F, M, S, p)
s = plotStyle(); nS = size(S.Jbest, 2);
fig = newFigure('fleetSweep', [50 50 1400 950]);
tl = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf(['Problem 4: fixed-phase sweep at \\tau = t_f (eq. (65)), joint optimum ' ...
    '\\phi_\\tau^* = %.2f\\circ, J_4^* = %.2f N s'], rad2deg(M.phiTau), M.J4));
g = ~S.isPresc;
[x, o] = sort(rad2deg(S.phis(g)));
Jg = S.Jbest(g, :); Jg = Jg(o, :);
cg = S.clsBest(g, :); cg = cg(o, :);
dg = S.dJbest(g, :); dg = dg(o, :);
if x(1) == -180, x = [x, 180]; Jg = [Jg; Jg(1,:)]; cg = [cg; cg(1,:)]; dg = [dg; dg(1,:)]; end
phiS = rad2deg(M.phiTau); phiP = rad2deg(S.phiPresc);

% (a) stacked contributions, combiner at the bottom
ax = nexttile(tl, [1 3]); hold(ax, 'on'); grid(ax, 'on');
ord = [5 1:4];
ha = area(ax, x, Jg(:, ord), 'LineStyle', 'none');
for k = 1:nS
    ha(k).FaceColor = s.sc(ord(k),:); ha(k).FaceAlpha = 0.8;
    ha(k).DisplayName = scLabel(p, ord(k));
end
plot(ax, x, sum(Jg, 2), 'k.-', 'LineWidth', 1.5, 'MarkerSize', 14, 'DisplayName', 'J_4(\phi_\tau) = \Sigma_i J_i(t_f, \phi_\tau)');
for k = 1:numel(F.starts)                % joint starts: grid point -> joint optimum
    st = F.starts(k);
    if ~strcmp(st.status, 'Solve_Succeeded'), continue, end
    plot(ax, [st.phiGridDeg, st.phiTauDeg], [st.J_grid_Ns, st.J_Ns], '.-', 'Color', [0.35 0.35 0.35], ...
         'MarkerSize', 10, 'LineWidth', 1, 'HandleVisibility', 'off');
end
plot(ax, nan, nan, '.-', 'Color', [0.35 0.35 0.35], 'DisplayName', ...
     'joint NLP: best grid points \rightarrow joint optima');
plot(ax, phiS, M.J4, 'p', 'MarkerSize', 18, 'MarkerFaceColor', 'y', 'Color', 'k', ...
     'DisplayName', sprintf('joint optimum: \\phi_\\tau^* = %.2f\\circ, J_4^* = %.2f N s', phiS, M.J4));
plot(ax, phiP, M.prescBound, 's', 'MarkerSize', 10, 'MarkerFaceColor', s.prob(2,:), 'Color', 'k', ...
     'DisplayName', sprintf('Problem 2 phase %.0f\\circ: %.2f N s (upper bound)', phiP, M.prescBound));
for i = 1:nS
    xline(ax, rad2deg(angle(exp(1i*M.phi3i(i)))), '-', 'Color', s.sc(i,:), 'LineWidth', 1.5, ...
          'HandleVisibility', 'off');
end
lo = min(sum(Jg, 2)); hi = max(sum(Jg, 2));
ylim(ax, [0, hi*1.08]); xlim(ax, [-180 180]); xticks(ax, -180:45:180);
xlabel(ax, 'Arrival phase \phi_\tau (deg)'); ylabel(ax, 'Impulse (N s)');
legend(ax, 'Location', 'eastoutside');
title(ax, sprintf(['(a) per-spacecraft costs J_i(t_f, \\phi_\\tau), stacked; vertical lines: individually ' ...
    'optimal phases \\phi^*_{3,i} (Problem 3); J_4 range %.0f-%.0f N s'], lo, hi), 'FontWeight', 'normal');

% (b) where each spacecraft pulls: J_i - min J_i, with the multiplier slopes
ax = nexttile(tl, [1 2]); hold(ax, 'on'); grid(ax, 'on');
for i = 1:nS
    y = Jg(:, i) - min(Jg(:, i));
    plot(ax, x, y, '-', 'Color', s.sc(i,:), 'LineWidth', 1.5, 'DisplayName', scLabel(p, i));
    n = cg(:, i) > 0;
    plot(ax, x(n), y(n), '^', 'Color', s.sc(i,:), 'MarkerFaceColor', s.sc(i,:), 'HandleVisibility', 'off');
    plot(ax, x(~n), y(~n), 'v', 'Color', s.sc(i,:), 'MarkerFaceColor', 'w', 'HandleVisibility', 'off');
    for k = 1:numel(x)                   % tangents from the phase multiplier
        dx = 6; sl = dg(k, i)*pi/180;
        plot(ax, x(k) + [-dx dx], y(k) + sl*[-dx dx], '-', 'Color', [s.sc(i,:) 0.6], ...
             'LineWidth', 0.8, 'HandleVisibility', 'off');
    end
end
xline(ax, phiS, 'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
xlim(ax, [-180 180]); xticks(ax, -180:45:180);
xlabel(ax, 'Arrival phase \phi_\tau (deg)'); ylabel(ax, 'J_i - min_\phi J_i (N s)');
legend(ax, 'Location', 'north', 'NumColumns', 3);
title(ax, '(b) per spacecraft; \^ north, v south path class; short lines: slope dJ_i/d\phi from the multipliers; dashed: \phi_\tau^*', ...
      'FontWeight', 'normal');

% (c) sensitivities at the optimum
ax = nexttile(tl); hold(ax, 'on'); grid(ax, 'on');
b = bar(ax, 1:nS, M.dJdphi, 'FaceColor', 'flat');
b.CData = s.sc;
plot(ax, 1:nS, M.dJdphiFD, 'kd', 'MarkerFaceColor', 'w', 'MarkerSize', 8);
yline(ax, 0, 'k-', 'HandleVisibility', 'off');
xticks(ax, 1:nS); xticklabels(ax, compose('%d', 1:nS)); xlabel(ax, 'Spacecraft');
ylabel(ax, '\partial J_i/\partial\phi_\tau at \phi_\tau^* (N s/rad)');
legend(ax, {'terminal multipliers', sprintf('finite differences (fit over \\pm%g\\circ)', 2*F.fd.stepDeg)}, ...
       'Location', 'best');
title(ax, sprintf('(c) sensitivities, sum %.1e N s/rad', M.sensSum), 'FontWeight', 'normal');
end

function fig = impulseFig(F, M, p)
s = plotStyle(); nS = numel(F.res);
fig = newFigure('fleetImpulse', [50 50 1400 600]);
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf(['Problem 4: J_4^* = %.2f N s (\\Deltam = %.3f kg), \\Sigma_i J^*_{3,i} = %.2f N s, ' ...
    'price of the shared \\tau, \\phi_\\tau: \\Pi_4 = %.2f N s (%.2f%%)'], M.J4, M.dm, sum(M.J3i), ...
    M.Pi4, 100*M.Pi4rel));
ax = nexttile(tl); hold(ax, 'on'); grid(ax, 'on');
b = bar(ax, 1:nS, [M.Ji; M.J3i]', 'grouped', 'FaceColor', 'flat');
b(1).CData = s.sc; b(2).CData = 0.45 + 0.55*s.sc;
for i = 1:nS
    text(ax, i, M.Ji(i), sprintf('%.1f%%\n%.1f g', 100*M.share(i), 1e3*M.dmi(i)), ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 9);
end
xticks(ax, 1:nS); xticklabels(ax, arrayfun(@(i) scLabel(p, i), 1:nS, 'UniformOutput', false));
ylabel(ax, 'Impulse (N s)'); ylim(ax, [0 1.18*max(M.Ji)]);
legend(ax, {'J_{4,i}: shared \tau^*, \phi_\tau^*', 'J^*_{3,i}: own \tau, \phi_\tau (Problem 3)'}, ...
       'Location', 'northwest');
title(ax, 'Impulse per spacecraft (labels: share of J_4^*, propellant)', 'FontWeight', 'normal');
ax = nexttile(tl); hold(ax, 'on'); grid(ax, 'on');
tot = 0;
for i = 1:nS
    r = F.res{i};
    plot(ax, r.sol.tNode, r.cumImpulse, '-', 'Color', s.sc(i,:), 'LineWidth', 1.6, ...
         'DisplayName', scLabel(p, i));
    tot = tot + r.cumImpulse;
end
plot(ax, F.res{1}.sol.tNode, tot, 'k--', 'LineWidth', 1.6, 'DisplayName', 'fleet');
xlim(ax, [0 M.tau]); xlabel(ax, 'Time t (s)'); ylabel(ax, 'Cumulative impulse (N s)');
legend(ax, 'Location', 'northwest');
title(ax, 'Cumulative impulse', 'FontWeight', 'normal');
end

function fig = trajectoryFig(F, p)
s = plotStyle(); nS = numel(F.res);
fig = newFigure('fleetTraj', [50 50 1500 850]);
tl = tiledlayout(fig, 2, 5, 'TileSpacing', 'compact', 'Padding', 'compact');
pk = F.res{1}.p;
title(tl, sprintf(['Problem 4 paths in the target formation frame O^f (\\phi_\\tau^* = %.2f\\circ): ' ...
    'target pyramid red, initial pyramid at t = 0 grey; arrows: z^i every 40 min'], ...
    rad2deg(F.res{1}.sol.phiTau)));
ax = nexttile(tl, [2 2]); hold(ax, 'on'); grid(ax, 'on');
drawFormation(ax, [pk.viz.collectors, pk.viz.combiner]);
drawFormation(ax, initialFormationPoints(pk), [0.55 0.55 0.55]);
h = gobjects(1, nS);
for i = 1:nS
    r = F.res{i}; rf = r.prop.x(8:10,:) + r.p.rfp;
    h(i) = plot3(ax, rf(1,:), rf(2,:), rf(3,:), '-', 'Color', s.sc(i,:), 'LineWidth', 1.8);
    plot3(ax, rf(1,1), rf(2,1), rf(3,1), 's', 'MarkerFaceColor', s.sc(i,:), 'Color', 'k', 'MarkerSize', 8);
    arrows(ax, r, 60, s.sc(i,:));
end
axis(ax, 'equal'); view(ax, 35, 20);
xlabel(ax, 'x_f (m)'); ylabel(ax, 'y_f (m)'); zlabel(ax, 'z_f (m)');
legend(ax, h, arrayfun(@(i) scLabel(p, i), 1:nS, 'UniformOutput', false), 'Location', 'northwest');
title(ax, '(a) all spacecraft (squares: start)', 'FontWeight', 'normal');
% close-ups: collectors seen along z_f on a common scale, combiner along y_f
ext = 0;
for i = 1:min(nS, 4)
    d = F.res{i}.prop.x(8:10,:);
    ext = max(ext, max(abs(d(:))));
end
tiles = [3 4 8 9 5];
dims = {[1 2], [1 2], [1 2], [1 2], [1 3]};
for i = 1:nS
    r = F.res{i}; rf = r.prop.x(8:10,:) + r.p.rfp; c = r.p.rfp; dd = dims{i};
    if i < 5, ax = nexttile(tl, tiles(i)); else, ax = nexttile(tl, tiles(i), [2 1]); end
    hold(ax, 'on'); grid(ax, 'on');
    plot(ax, rf(dd(1),:), rf(dd(2),:), '-', 'Color', s.sc(i,:), 'LineWidth', 1.6);
    plot(ax, rf(dd(1),1), rf(dd(2),1), 's', 'MarkerFaceColor', s.sc(i,:), 'Color', 'k');
    plot(ax, c(dd(1)), c(dd(2)), 'o', 'MarkerSize', 9, 'MarkerFaceColor', [0.8 0.1 0.1], 'Color', 'k');
    if i < 5
        a = p.target.nulling/2; b = p.target.imaging/2;      % rectangle edges
        plot(ax, [a -a -a a a], [b b -b -b b], ':', 'Color', [0.8 0.1 0.1], 'LineWidth', 1.2);
        xlabel(ax, 'x_f (m)'); ylabel(ax, 'y_f (m)');
    else
        plot(ax, [0 0], [0 p.target.height], '--', 'Color', [0.8 0.1 0.1]);
        xlabel(ax, 'x_f (m)'); ylabel(ax, 'z_f (m)');
    end
    axis(ax, 'equal');
    if i < 5
        L = 1.05*ext; axis(ax, [c(1) - L, c(1) + L, c(2) - L, c(2) + L]);
    else
        pad = 60; axis(ax, [min(rf(1,:)) - pad, max(rf(1,:)) + pad, min(rf(3,:)) - pad, max(rf(3,:)) + pad]);
    end
    title(ax, sprintf('(%c) %s', 'b' + i - 1, scLabel(p, i)), 'FontWeight', 'normal', 'Color', s.sc(i,:));
end
end

function arrows(ax, r, len, col)
for t = 0:2400:r.sol.tau
    xs = stateAt(r, t);
    R = quatToDcm(xs(1:4)); z = R(:,3); o = xs(8:10) + r.p.rfp;
    quiver3(ax, o(1), o(2), o(3), len*z(1), len*z(2), len*z(3), 0, 'Color', col, ...
            'LineWidth', 1.2, 'MaxHeadSize', 0.6, 'HandleVisibility', 'off');
end
end

function fig = trackingFig(F, p)
s = plotStyle(); nS = numel(F.res);
fig = newFigure('fleetTracking');
tau = F.res{1}.sol.tau;
tl = tiledlayout(fig, 4, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, 'Problem 4 target-relative errors (independent propagation), norms per spacecraft');
spec = {@(r) r.thetaAcos, '\xi_i (rad)', p.gate.theta
        @(r) vecnorm(r.prop.x(5:7,:)), '|\omega^f_{fi}| (1/s)', p.gate.w
        @(r) vecnorm(r.prop.x(8:10,:)), '|r^f_{pi}| (m)', p.gate.r
        @(r) vecnorm(r.prop.x(11:13,:)), '|dr^f_{pi}/dt| (m/s)', p.gate.v};
t = F.res{1}.prop.t; zoom = t >= 0.9*tau;
for k = 1:4
    ax = nexttile(tl, [1 2]); hold(ax, 'on'); grid(ax, 'on');
    ax2 = nexttile(tl); hold(ax2, 'on'); grid(ax2, 'on');
    lo = inf;
    for i = 1:nS
        y = max(spec{k, 1}(F.res{i}), 1e-14);
        plot(ax, t, y, '-', 'Color', s.sc(i,:), 'LineWidth', 1.2);
        plot(ax2, t(zoom), y(zoom), '-', 'Color', s.sc(i,:), 'LineWidth', 1.2);
        plot(ax2, tau, y(end), 'o', 'MarkerFaceColor', s.sc(i,:), 'Color', 'k');
        lo = min(lo, y(end));
    end
    yline(ax2, spec{k, 3}, 'g-', 'tolerance', 'LineWidth', 1.2, 'LabelHorizontalAlignment', 'left');
    set([ax ax2], 'YScale', 'log');
    ylabel(ax, spec{k, 2}); xlim(ax, [0 tau]); xlim(ax2, [t(find(zoom, 1)) tau]);
    ylim(ax2, [min(lo, spec{k, 3})/10, inf]);
    if k == 1
        legend(ax, arrayfun(@(i) scLabel(p, i), 1:nS, 'UniformOutput', false), ...
               'Location', 'southwest', 'NumColumns', 5);
        title(ax2, 'arrival (last 10%) and gate', 'FontWeight', 'normal');
    end
    if k == 4, xlabel(ax, 'Time t (s)'); xlabel(ax2, 'Time t (s)'); end
end
end

function fig = sunFig(F, M, p)
s = plotStyle(); nS = numel(F.res);
fig = newFigure('fleetSun', [50 50 1500 650]);
tl = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, 'Problem 4 Sun-angle constraint (48) along the propagated trajectories');
ax = nexttile(tl, [1 2]); hold(ax, 'on'); grid(ax, 'on');
tau = M.tau;
sunBands(ax, p, [0 tau], struct('t0', {}, 't1', {}), s);
h = gobjects(1, nS); cl = ["south", "", "north"];
for i = 1:nS
    sn = F.res{i}.sun;
    h(i) = plot(ax, sn.t, sn.antiDeg, '-', 'Color', s.sc(i,:), 'LineWidth', 1.6, ...
        'DisplayName', sprintf('%s: %s, max violation %.4f\\circ, active at %d points', ...
        scLabel(p, i), cl(M.cls(i) + 2), sn.maxViolDeg, numel(sn.active)));
    if ~isempty(sn.active)
        plot(ax, [sn.active.t], [sn.active.antiDeg], 'o', 'MarkerFaceColor', s.sc(i,:), ...
             'Color', 'k', 'HandleVisibility', 'off');
    end
end
xlim(ax, [0 tau]); ylim(ax, [30 100]); yticks(ax, 30:10:100);
xlabel(ax, 'Time t (s)'); ylabel(ax, 'Anti-Sun angle of z^i (deg)');
legend(ax, h, 'Location', 'southoutside');
title(ax, 'green: field of regard; dots: constraint active', 'FontWeight', 'normal');
ax = nexttile(tl); drawSky(ax, p);
for i = 1:nS
    P = skyXY(F.res{i}.sun);
    plot(ax, P(1,:), P(2,:), '-', 'Color', s.sc(i,:), 'LineWidth', 1.6);
    plot(ax, P(1,1), P(2,1), 's', 'MarkerFaceColor', s.sc(i,:), 'Color', 'k');
    plot(ax, P(1,end), P(2,end), 'o', 'MarkerFaceColor', s.sc(i,:), 'Color', 'k');
end
title(ax, 'z^i about the Sun line (square: start, circle: arrival)', 'FontWeight', 'normal');
end

function fig = thrustFig(F, p)
s = plotStyle(); nS = numel(F.res);
fig = newFigure('fleetThrust', [50 50 1400 1000]);
tl = tiledlayout(fig, nS, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
uMax = max(cellfun(@(r) max(r.sol.U(:)), F.res));
title(tl, sprintf(['Problem 4 firing sequences T_{ik}/T_{max} (rows: thrusters by pod; colour scale to ' ...
    'the fleet peak %.2f) and total thrust \\Sigma_k T_{ik}'], uMax));
for i = 1:nS
    sol = F.res{i}.sol;
    ax = nexttile(tl);
    C = [sol.U; sol.U(end,:)]; C = [C, C(:,end)];
    [X, Y] = meshgrid(sol.tNode, 0.5:1:p.nThr + 0.5);
    surface(ax, X, Y, zeros(size(X)), C, 'EdgeColor', 'none');
    colormap(ax, flipud(hot)); clim(ax, [0 uMax]);
    hold(ax, 'on');
    for b = 5.5:5:15.5, yline(ax, b, '-', 'Color', [0.6 0.6 0.6]); end
    set(ax, 'YDir', 'reverse', 'YTick', 3:5:18, 'YTickLabel', compose('pod %d', 1:4), 'FontSize', 8);
    ylim(ax, [0.5 p.nThr + 0.5]); xlim(ax, [0 sol.tau]); box(ax, 'on');
    yyaxis(ax, 'right');
    stairs(ax, sol.tNode, [sum(sol.T, 1), sum(sol.T(:, end))], '-', 'Color', s.sc(i,:), 'LineWidth', 1.2);
    ylabel(ax, '\Sigma_k T (N)'); ax.YAxis(2).Color = s.sc(i,:);
    yyaxis(ax, 'left');
    ylabel(ax, sprintf('%d', i), 'Color', s.sc(i,:), 'FontWeight', 'bold', 'FontSize', 12);
    title(ax, sprintf('%s: J = %.2f N s, peak %.2f N', scLabel(p, i), F.res{i}.J, F.res{i}.maxThrust), ...
          'FontWeight', 'normal', 'Color', s.sc(i,:));
    if i == nS, xlabel(ax, 'Time t (s)'); else, xticklabels(ax, {}); end
end
cb = colorbar(ax); cb.Layout.Tile = 'east'; cb.Label.String = 'T/T_{max}';
end

function fig = distanceFig(M, ~)
fig = newFigure('fleetDistance', [50 50 1300 600]);
ax = axes(fig); hold(ax, 'on'); grid(ax, 'on');
d = M.dist;
cmb = any(d.pairs == 5, 2);
for k = 1:size(d.pairs, 1)
    col = [0.6 0.6 0.6]; if cmb(k), col = [0.75 0.6 0.8]; end
    plot(ax, d.t, d.D(k,:), '-', 'Color', col, 'LineWidth', 1, 'HandleVisibility', 'off');
    text(ax, d.t(end), d.D(k, end), sprintf(' %d-%d', d.pairs(k,:)), 'FontSize', 7, 'Color', 0.6*col);
end
plot(ax, d.t, min(d.D, [], 1), 'k-', 'LineWidth', 2.2, 'DisplayName', 'minimum over all pairs');
plot(ax, d.minTime, d.min, 'rp', 'MarkerSize', 14, 'MarkerFaceColor', 'r', 'DisplayName', ...
     sprintf('closest approach %.1f m: spacecraft %d and %d at t = %.0f s', d.min, d.minPair, d.minTime));
set(ax, 'YScale', 'log'); xlim(ax, [0 d.t(end)]);
xlabel(ax, 'Time t (s)'); ylabel(ax, 'Distance between spacecraft centres (m)');
legend(ax, 'Location', 'southwest');
title(ax, sprintf(['Problem 4: distances between all %d pairs (grey: collectors, purple: with the ' ...
    'combiner); not constrained (Assumption 5), reference for collision avoidance'], ...
    size(d.pairs, 1)), 'FontWeight', 'normal');
end
