function fig = plotThrusters(res, p)
%PLOTTHRUSTERS Twenty thrust magnitudes as step histories, one panel per
%   pod, with the bounds 0 and Tmax (nonuniform intervals to scale).
sol = res.sol;
fig = newFigure('thrusters');
tl = tiledlayout(fig, 4, 1, 'TileSpacing', 'compact');
title(tl, 'Thrust magnitudes T_{ik} grouped by pod');
cols = lines(5);
for j = 1:4
    ax = nexttile(tl); hold(ax, 'on'); grid(ax, 'on');
    ks = find(p.thr.pod == j);
    for c = 1:numel(ks)
        k = ks(c);
        stairs(ax, sol.tNode, [sol.T(k,:), sol.T(k,end)], 'Color', cols(c,:), 'LineWidth', 1.2);
    end
    yline(ax, [0 p.Tmax], 'k--');
    ylim(ax, [-0.5 p.Tmax + 0.5]); xlim(ax, [0 sol.tau]);
    ylabel(ax, sprintf('Pod %d (N)', j));
    legend(ax, cellstr(p.thr.label(ks)), 'Location', 'eastoutside');
end
xlabel(ax, 'Time t (s)');
end
