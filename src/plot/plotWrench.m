function fig = plotWrench(res, ~)
%PLOTWRENCH Body-frame resultant control force and torque, (29)-(30).
s = plotStyle(); sol = res.sol;
fig = newFigure('wrench', [50 50 1300 650]);
tl = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact');
title(tl, 'Resultant control wrench in the body frame O^i');
Y = {res.f, 'f^i_{ctrl} (N)'; res.tau, '\tau^i_{ctrl} (N m)'};
for i = 1:2
    ax = nexttile(tl); hold(ax, 'on'); grid(ax, 'on');
    for c = 1:3
        stairs(ax, sol.tNode, [Y{i,1}(c,:), Y{i,1}(c,end)], 'Color', s.comp(c,:), 'LineWidth', 1.2);
    end
    ylabel(ax, Y{i,2}); xlim(ax, [0 sol.tau]);
    legend(ax, {'x^i', 'y^i', 'z^i'}, 'Location', 'northeast');
end
xlabel(ax, 'Time t (s)');
end
