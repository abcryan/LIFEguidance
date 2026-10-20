function fig = plotHeatmap(res, p)
%PLOTHEATMAP T_ik/Tmax over (nonuniform) intervals, full span and burns.
sol = res.sol;
fig = newFigure('heatmap', [50 50 1300 600]);
burns = res.burns;
tl = tiledlayout(fig, 2, max(numel(burns), 1), 'TileSpacing', 'compact');
title(tl, 'Firing sequence T_{ik}/T_{max}');
ax = nexttile(tl, [1 max(numel(burns), 1)]);
drawMap(ax, sol, p, [0 sol.tau]);
cb = colorbar(ax); cb.Label.String = 'T/T_{max}';
for b = 1:numel(burns)
    ax = nexttile(tl);
    pad = max(0.25*burns(b).duration, 2);
    drawMap(ax, sol, p, [burns(b).t0 - pad, burns(b).t1 + pad]);
    title(ax, sprintf('burn %d', b), 'FontWeight', 'normal');
end
end

function drawMap(ax, sol, p, xl)
C = [sol.U; sol.U(end,:)]; C = [C, C(:,end)];
[X, Y] = meshgrid(sol.tNode, 0.5:1:p.nThr + 0.5);
surface(ax, X, Y, zeros(size(X)), C, 'EdgeColor', 'none');
colormap(ax, flipud(hot)); clim(ax, [0 1]);
set(ax, 'YDir', 'reverse', 'YTick', 1:p.nThr, 'YTickLabel', p.thr.label, 'FontSize', 7);
xlim(ax, xl); ylim(ax, [0.5 p.nThr + 0.5]); xlabel(ax, 'Time t (s)'); box(ax, 'on');
end
