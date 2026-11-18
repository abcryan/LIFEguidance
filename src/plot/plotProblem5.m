function figs = plotProblem5(F, M, P4, p, outDir)
%PLOTPROBLEM5 Problem 5 figures (Sec. 5.3), saved as PNG to outDir:
%   p5_separation   d_ij(t) of the crossing pairs (Problem 5 solid,
%                   Problem 4 dashed) with dmin, the NLP margin and the
%                   active intervals; all ten pairs on a log scale
%   p5_passes       relative trajectory of each crossing pair in the frame
%                   of its first spacecraft (axes parallel to O*), with the
%                   keep-out sphere of radius dmin and the closest approaches
%   p5_collectors   collector paths in O* (Problem 4 grey), with the
%                   spheres of radius dmin/2 at the closest approaches
%   p5_impulse      impulse per spacecraft, Problems 4 and 5, Pi5
%   p5_total_cost   total fleet impulse of Problems 4 and 5 with Pi5, and the
%                   change per spacecraft
%   p5_dmin_sweep   J5(dmin) with J4*, Problem 4's minimum separation and
%                   the multiplier tangent at dmin
list = {
    'p5_separation', @() separationFig(F, M, p)
    'p5_passes',     @() passesFig(F, M, P4, p)
    'p5_collectors', @() collectorsFig(F, M, P4, p)
    'p5_impulse',    @() impulseFig(F, M, p)
    'p5_total_cost', @() totalCostFig(F, M, P4, p)
    'p5_dmin_sweep', @() sweepFig(F, M, p)
    };
for k = 1:size(list, 1)
    figs.(list{k, 1}) = list{k, 2}();
    exportgraphics(figs.(list{k, 1}), fullfile(outDir, [list{k, 1} '.png']), 'Resolution', 200);
end
end

function cp = crossPairs(F)
cp = reshape([F.crossing.pair], 2, [])';
end

function fig = separationFig(F, M, p)
s = plotStyle(); cd = p.collision;
fig = newFigure('p5Separation', [50 50 1400 900]);
tl = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf(['Problem 5: separations d_{ij}(t) (propagated); collision avoidance d_{min} = %.0f m, ' ...
    'imposed with %.1f m in the NLP'], cd.dmin, cd.dmin + cd.margin));
ax = nexttile(tl); hold(ax, 'on'); grid(ax, 'on');
for a = M.activeIntervals
    patch(ax, [a.t0 a.t1 a.t1 a.t0], [0 0 80 80], [1 0.85 0.85], 'EdgeColor', 'none', ...
          'HandleVisibility', 'off');
end
d = M.dist; d4 = M.dist4; cps = crossPairs(F);
h = gobjects(0);
for k = 1:size(cps, 1)
    kk = find(ismember(d.pairs, cps(k,:), 'rows')); c = s.sc(cps(k,1), :);
    h(end+1) = plot(ax, d.t, d.D(kk,:), '-', 'Color', c, 'LineWidth', 2, 'DisplayName', ...
        sprintf('(%d,%d) Problem 5: min %.2f m', cps(k,:), M.minima(kk).dMin)); %#ok<AGROW>
    [m4, i4] = min(d4.D(kk,:));
    h(end+1) = plot(ax, d4.t, d4.D(kk,:), '--', 'Color', c, 'LineWidth', 1.4, 'DisplayName', ...
        sprintf('(%d,%d) Problem 4: min %.2f m', cps(k,:), m4)); %#ok<AGROW>
    plot(ax, d4.t(i4), m4, 'x', 'Color', c, 'MarkerSize', 12, 'LineWidth', 2, 'HandleVisibility', 'off');
    plot(ax, M.minima(kk).tMin, M.minima(kk).dMin, 'o', 'MarkerFaceColor', c, 'Color', 'k', ...
         'MarkerSize', 8, 'HandleVisibility', 'off');
end
h(end+1) = yline(ax, cd.dmin, 'r-', 'LineWidth', 1.6, 'DisplayName', sprintf('d_{min} = %.0f m', cd.dmin));
h(end+1) = yline(ax, cd.dmin + cd.margin, 'r:', 'LineWidth', 1.2, 'DisplayName', 'NLP bound d_{min} + 0.5 m');
h(end+1) = patch(ax, nan, nan, [1 0.85 0.85], 'EdgeColor', 'none', 'DisplayName', 'constraint active (collocation)');
xlim(ax, [0 d.t(end)]); ylim(ax, [0 80]);
ylabel(ax, 'd_{ij} (m)'); legend(ax, h, 'Location', 'northwest', 'NumColumns', 2);
title(ax, 'crossing pairs (collectors that swap sides); x: Problem 4 closest approach, o: Problem 5', ...
      'FontWeight', 'normal');
ax = nexttile(tl); hold(ax, 'on'); grid(ax, 'on');
for k = 1:size(d.pairs, 1)
    isC = ismember(d.pairs(k,:), cps, 'rows');
    col = [0.6 0.6 0.6]; if any(d.pairs(k,:) == 5), col = [0.75 0.6 0.8]; end
    if isC, col = s.sc(d.pairs(k,1), :); end
    plot(ax, d.t, d.D(k,:), '-', 'Color', col, 'LineWidth', 1 + isC, 'HandleVisibility', 'off');
    plot(ax, d4.t, d4.D(k,:), ':', 'Color', col, 'LineWidth', 0.8, 'HandleVisibility', 'off');
    text(ax, d.t(end), d.D(k, end), sprintf(' %d-%d', d.pairs(k,:)), 'FontSize', 7, 'Color', 0.6*col);
end
yline(ax, cd.dmin, 'r-', 'LineWidth', 1.4);
set(ax, 'YScale', 'log'); xlim(ax, [0 d.t(end)]);
xlabel(ax, 'Time t (s)'); ylabel(ax, 'd_{ij} (m)');
title(ax, sprintf(['all %d pairs (solid: Problem 5, dotted: Problem 4; grey: collectors, purple: ' ...
    'with the combiner); smallest Problem 5 separation %.3f m (%d,%d)'], size(d.pairs, 1), ...
    M.minSep, M.minSepPair), 'FontWeight', 'normal');
end

function fig = passesFig(F, ~, P4, p)
s = plotStyle(); cd = p.collision; cps = crossPairs(F); nP = size(cps, 1);
fig = newFigure('p5Passes', [50 50 1400 650*nP]);
tl = tiledlayout(fig, nP, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf(['Problem 5: relative motion of each crossing pair, r_j - r_i in the frame of ' ...
    'spacecraft i (axes parallel to O^*); red sphere: keep-out radius d_{min} = %.0f m'], cd.dmin));
[xs, ys, zs] = sphere(40);
for k = 1:nP
    i = cps(k,1); j = cps(k,2);
    r5 = F.res{j}.prop.xP(8:10,:) - F.res{i}.prop.xP(8:10,:);
    r4 = P4.fleet.res{j}.prop.xP(8:10,:) - P4.fleet.res{i}.prop.xP(8:10,:);
    v5 = F.res{j}.prop.xP(11:13,:) - F.res{i}.prop.xP(11:13,:);
    [d5, i5] = min(vecnorm(r5)); [d4, i4] = min(vecnorm(r4));
    n = cross(r5(:, i5), v5(:, i5)); n = n/norm(n);   % plane of the closest approach
    for v = 1:2
        ax = nexttile(tl); hold(ax, 'on'); grid(ax, 'on'); axis(ax, 'equal');
        surf(ax, cd.dmin*xs, cd.dmin*ys, cd.dmin*zs, 'FaceColor', [1 0.3 0.3], 'FaceAlpha', 0.15, ...
             'EdgeColor', [1 0.5 0.5], 'EdgeAlpha', 0.25, 'DisplayName', sprintf('keep-out, %.0f m', cd.dmin));
        plot3(ax, 0, 0, 0, 'k+', 'MarkerSize', 10, 'LineWidth', 1.5, 'HandleVisibility', 'off');
        plot3(ax, r4(1,:), r4(2,:), r4(3,:), '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 1.4, ...
              'DisplayName', sprintf('Problem 4: closest %.2f m', d4));
        plot3(ax, r5(1,:), r5(2,:), r5(3,:), '-', 'Color', s.sc(j,:), 'LineWidth', 2, ...
              'DisplayName', sprintf('Problem 5: closest %.2f m', d5));
        plot3(ax, r4(1,i4), r4(2,i4), r4(3,i4), 'x', 'Color', 'k', 'MarkerSize', 12, 'LineWidth', 2, ...
              'HandleVisibility', 'off');
        plot3(ax, r5(1,i5), r5(2,i5), r5(3,i5), 'o', 'MarkerFaceColor', s.sc(j,:), 'Color', 'k', ...
              'MarkerSize', 8, 'HandleVisibility', 'off');
        plot3(ax, [0 r5(1,i5)], [0 r5(2,i5)], [0 r5(3,i5)], '-', 'Color', s.sc(j,:), 'HandleVisibility', 'off');
        plot3(ax, r5(1,1), r5(2,1), r5(3,1), 's', 'MarkerFaceColor', 'w', 'Color', 'k', 'HandleVisibility', 'off');
        plot3(ax, r5(1,end), r5(2,end), r5(3,end), 'p', 'MarkerFaceColor', 'w', 'Color', 'k', ...
              'MarkerSize', 12, 'HandleVisibility', 'off');
        L = 1.2*max(vecnorm([r4(:, [1 end]), r5(:, [1 end])]));
        lim = max(L, 3*cd.dmin);
        xlim(ax, [-lim lim]); ylim(ax, [-lim lim]); zlim(ax, [-lim lim]);
        xlabel(ax, 'x_* (m)'); ylabel(ax, 'y_* (m)'); zlabel(ax, 'z_* (m)');
        if v == 1
            view(ax, 35, 20);
            title(ax, sprintf('(%d,%d): spacecraft %d seen from %d (square: t = 0, star: arrival)', ...
                  i, j, j, i), 'FontWeight', 'normal');
            legend(ax, 'Location', 'northeast');
        else
            view(ax, n'); lim2 = 3*cd.dmin;
            xlim(ax, [-lim2 lim2]); ylim(ax, [-lim2 lim2]); zlim(ax, [-lim2 lim2]);
            title(ax, 'zoom, viewed normal to the plane of the Problem 5 closest approach', ...
                  'FontWeight', 'normal');
        end
        camlight(ax, 'headlight'); lighting(ax, 'gouraud');
    end
end
end

function fig = collectorsFig(F, M, P4, p)
s = plotStyle(); cd = p.collision; cps = crossPairs(F);
fig = newFigure('p5Collectors', [50 50 1500 800]);
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf(['Problem 5: collector paths in the pointing frame O^* (grey: Problem 4); circles of ' ...
    'radius d_{min}/2 = %.0f m at the closest approach of each crossing pair touch exactly at d_{min}'], cd.dmin/2));
views = {[1 2], [1 3]}; lab = {'x_* (m)', 'y_* (m)', 'z_* (m)'};
ph = linspace(0, 2*pi, 90);
for v = 1:2
    ax = nexttile(tl); hold(ax, 'on'); grid(ax, 'on'); axis(ax, 'equal');
    a = views{v}(1); b = views{v}(2);
    for i = 1:4
        r4 = P4.fleet.res{i}.prop.xP(8:10,:); r5 = F.res{i}.prop.xP(8:10,:);
        plot(ax, r4(a,:), r4(b,:), '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2, 'HandleVisibility', 'off');
        plot(ax, r5(a,:), r5(b,:), '-', 'Color', s.sc(i,:), 'LineWidth', 1.8, 'DisplayName', ...
             sprintf('%d %s', i, p.fleet.name(i)));
        plot(ax, r5(a,1), r5(b,1), 's', 'MarkerFaceColor', s.sc(i,:), 'Color', 'k', 'HandleVisibility', 'off');
        plot(ax, r5(a,end), r5(b,end), 'p', 'MarkerFaceColor', s.sc(i,:), 'Color', 'k', 'MarkerSize', 11, ...
             'HandleVisibility', 'off');
    end
    for k = 1:size(cps, 1)
        kk = find(ismember(M.dist.pairs, cps(k,:), 'rows'));
        [~, i5] = min(M.dist.D(kk,:)); [~, i4] = min(M.dist4.D(kk,:));
        for e = cps(k,:)
            q4 = P4.fleet.res{e}.prop.xP(8:10, i4); q5 = F.res{e}.prop.xP(8:10, i5);
            plot(ax, q4(a) + cd.dmin/2*cos(ph), q4(b) + cd.dmin/2*sin(ph), '-', 'Color', [0.55 0.55 0.55], ...
                 'HandleVisibility', 'off');
            fill(ax, q5(a) + cd.dmin/2*cos(ph), q5(b) + cd.dmin/2*sin(ph), s.sc(e,:), 'FaceAlpha', 0.25, ...
                 'EdgeColor', s.sc(e,:), 'HandleVisibility', 'off');
        end
    end
    xlabel(ax, lab{a}); ylabel(ax, lab{b});
    if v == 1
        legend(ax, 'Location', 'eastoutside');
        title(ax, 'top view x_*-y_* (square: t = 0, star: arrival)', 'FontWeight', 'normal');
    else
        title(ax, 'side view x_*-z_* (grey circles: Problem 4 at its closest approach)', 'FontWeight', 'normal');
    end
end
end

function fig = impulseFig(F, M, p)
s = plotStyle(); nS = numel(F.res);
fig = newFigure('p5Impulse', [50 50 1300 550]);
ax = axes(fig); hold(ax, 'on'); grid(ax, 'on');
b = bar(ax, 1:nS, [M.J4i; M.Ji]', 'grouped');
b(1).FaceColor = [0.7 0.7 0.7]; b(2).FaceColor = 'flat'; b(2).CData = s.sc(1:nS,:);
for i = 1:nS
    text(ax, i + 0.15, M.Ji(i), sprintf('%+.2f', M.Ji(i) - M.J4i(i)), 'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'bottom', 'FontSize', 9);
end
xticks(ax, 1:nS); xticklabels(ax, compose('%d %s', (1:nS)', p.fleet.name(:)));
ylabel(ax, 'impulse (N s)'); legend(ax, {'Problem 4', 'Problem 5'}, 'Location', 'northwest');
title(ax, sprintf(['Problem 5: J_5^* = %.2f N s, J_4^* = %.2f N s, \\Pi_5 = %.2f N s (%.2f%%); ' ...
    '\\Delta\\tau^* = %+.1f s, \\Delta\\phi_\\tau^* = %+.2f\\circ'], M.J5, M.J4, M.Pi5, 100*M.Pi5rel, ...
    M.dTau, M.dPhiDeg), 'FontWeight', 'normal');
end

function fig = totalCostFig(F, M, ~, p)
s = plotStyle(); nS = numel(F.res);
fig = newFigure('p5TotalCost', [50 50 1250 560]);
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf(['Total fleet impulse: Problem 4 (no collision avoidance) vs Problem 5 (d_{min} = %.0f m): ' ...
    '\\Pi_5 = J_5^* - J_4^* = %.2f N s (%.2f%%)'], p.collision.dmin, M.Pi5, 100*M.Pi5rel), 'FontSize', 13);
% (a) totals, zoomed axis
ax = nexttile(tl); hold(ax, 'on'); grid(ax, 'on');
tot = [M.J4, M.J5];
bar(ax, 1, tot(1), 0.5, 'FaceColor', [0.6 0.6 0.6]);
bar(ax, 2, tot(2), 0.5, 'FaceColor', [0.2 0.2 0.2]);
for k = 1:2
    text(ax, k, tot(k), sprintf('%.2f N s', tot(k)), 'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'bottom', 'FontWeight', 'bold', 'FontSize', 11);
end
J4best = min([M.J4, F.sweep([F.sweep.success]).J_Ns, F.lowFromP4.J_Ns]);
yl = [floor(min(J4best, tot(1)) - 8), ceil(tot(2) + 8)];
plot(ax, [1.25 2.6], tot([1 1]), 'r:', 'LineWidth', 1.2);
plot(ax, [2.45 2.6 2.6 2.45], [tot(1) tot(1) tot(2) tot(2)], 'r-', 'LineWidth', 1.6);
text(ax, 2.63, mean(tot), sprintf('\\Pi_5\n%+.2f N s\n(%.2f%%)', M.Pi5, 100*M.Pi5rel), 'Color', 'r', ...
     'FontWeight', 'bold', 'VerticalAlignment', 'middle');
if J4best < M.J4 - M.tolJ
    plot(ax, [0.6 1.4], [J4best J4best], '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.2);
    text(ax, 0.62, J4best, sprintf('better Problem 4 optimum\nfound later: %.2f N s', J4best), ...
         'FontSize', 8, 'VerticalAlignment', 'top', 'Color', [0.3 0.3 0.3]);
end
xticks(ax, 1:2); xticklabels(ax, {'Problem 4', 'Problem 5'}); xlim(ax, [0.4 3.1]); ylim(ax, yl);
ylabel(ax, 'total fleet impulse (N s)');
title(ax, sprintf('(a) total (y-axis from %d N s)', yl(1)), 'FontWeight', 'normal');
% (b) change per spacecraft
ax = nexttile(tl); hold(ax, 'on'); grid(ax, 'on');
dJ = M.Ji - M.J4i; dJ(abs(dJ) < 5e-3) = 0;     % no "-0.00" 
for i = 1:nS
    bar(ax, i, dJ(i), 0.6, 'FaceColor', s.sc(i,:));
    text(ax, i, dJ(i), sprintf('%+.2f', dJ(i)), 'HorizontalAlignment', 'center', ...
         'VerticalAlignment', ifelse(dJ(i) >= 0, 'bottom', 'top'), 'FontSize', 10);
end
bar(ax, nS + 1.5, M.Pi5, 0.6, 'FaceColor', [0.2 0.2 0.2]);
text(ax, nS + 1.5, M.Pi5, sprintf('%+.2f', M.Pi5), 'HorizontalAlignment', 'center', ...
     'VerticalAlignment', 'bottom', 'FontWeight', 'bold', 'FontSize', 10);
yline(ax, 0, 'k-');
xticks(ax, [1:nS, nS + 1.5]); xticklabels(ax, [compose('%d', 1:nS), {'\Pi_5'}]);
xlim(ax, [0.3 nS + 2.2]); ylim(ax, [min(0, min(dJ)) - 2, M.Pi5 + 3]);
xlabel(ax, 'spacecraft (1-4 collectors, 5 combiner)'); ylabel(ax, 'J_{5,i} - J_{4,i} (N s)');
title(ax, '(b) who pays: change per spacecraft (sums to \Pi_5)', 'FontWeight', 'normal');
end

function s = ifelse(c, a, b)
if c, s = a; else, s = b; end
end

function fig = sweepFig(F, M, p)
cd = p.collision;
fig = newFigure('p5Sweep', [50 50 1100 650]);
ax = axes(fig); hold(ax, 'on'); grid(ax, 'on');
sw = F.sweep; ok = [sw.success];
plot(ax, [sw(ok).dmin], [sw(ok).J_Ns], 'o-', 'Color', [0.2 0.2 0.6], 'LineWidth', 1.6, ...
     'MarkerFaceColor', [0.2 0.2 0.6], 'DisplayName', 'J_5(d_{min})');
for k = find(ok)
    text(ax, sw(k).dmin, sw(k).J_Ns, sprintf('  %.2f N s/m', sw(k).dJddmin), 'FontSize', 8, ...
         'VerticalAlignment', 'top');
end
yline(ax, M.J4, 'k--', 'LineWidth', 1.2, 'DisplayName', sprintf('J_4^* = %.2f N s', M.J4));
xline(ax, M.dist4.min, ':', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.4, 'DisplayName', ...
      sprintf('Problem 4 min. separation %.2f m', M.dist4.min));
i0 = find([sw.dmin] == cd.dmin, 1);
dd = cd.dmin + [-4 4];
plot(ax, dd, sw(i0).J_Ns + M.dJddmin*(dd - cd.dmin), 'r-', 'LineWidth', 1.6, 'DisplayName', ...
     sprintf('multiplier tangent: %.3f N s/m', M.dJddmin));
plot(ax, F.fd.dmin, F.fd.J, 'rs', 'MarkerFaceColor', 'w', 'DisplayName', ...
     sprintf('central differences: %.3f N s/m', M.dJddminFD));
Jok = [sw(ok).J_Ns, M.J4]; yl = [min(Jok) - 5, max(Jok) + 8];
if any(~ok)                             % not converged: marked at the top, not at their value
    plot(ax, [sw(~ok).dmin], yl(2) - 2 + 0*[sw(~ok).dmin], 'x', 'Color', 'r', 'MarkerSize', 12, ...
         'LineWidth', 2, 'DisplayName', 'not converged');
    for k = find(~ok)
        text(ax, sw(k).dmin, yl(2) - 3.5, strrep(sw(k).status, '_', ' '), 'FontSize', 7, ...
             'HorizontalAlignment', 'center', 'Color', 'r');
    end
end
xline(ax, M.feasibility.initialMin - cd.margin, '-.', 'Color', [0.8 0.3 0.3], 'DisplayName', ...
      sprintf('upper limit: initial separation %.0f m - margin', M.feasibility.initialMin));
ylim(ax, yl); xlim(ax, [0 max([sw.dmin]) + 2]);
xlabel(ax, 'd_{min} (m)'); ylabel(ax, 'J_5^* (N s)'); legend(ax, 'Location', 'northwest');
title(ax, sprintf(['Problem 5: total impulse against the collision distance (continued from the %.0f m ' ...
    'solution; labels: dJ_5/dd_{min} from the multipliers)'], cd.dmin), 'FontWeight', 'normal');
end
