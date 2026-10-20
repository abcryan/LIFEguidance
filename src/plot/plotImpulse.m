function fig = plotImpulse(res, p)
%PLOTIMPULSE Cumulative total impulse and per-thruster integrated impulse.
sol = res.sol; s = plotStyle();
fig = newFigure('impulse', [50 50 1300 550]);
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact');
ax = nexttile(tl); grid(ax, 'on'); hold(ax, 'on');
plot(ax, sol.tNode, res.cumImpulse, 'k', 'LineWidth', 1.5);
xlabel(ax, 'Time t (s)'); ylabel(ax, 'Cumulative impulse (N s)');
title(ax, sprintf('Total %.3f N s  (\\Deltam = %.3f kg, I_{sp} = %g s)', res.J, res.dm, p.Isp), ...
      'FontWeight', 'normal');
xlim(ax, [0 sol.tau]);
ax = nexttile(tl); hold(ax, 'on'); grid(ax, 'on');
bh = bar(ax, 1:p.nThr, res.thrusterImpulse, 'FaceColor', 'flat');
bh.CData = s.pod(p.thr.pod, :);
set(ax, 'XTick', 1:p.nThr, 'XTickLabel', p.thr.label, 'XTickLabelRotation', 60);
ylabel(ax, 'Integrated impulse (N s)');
for k = 1:p.nThr
    if res.thrusterImpulse(k) > 0.05
        text(ax, k, res.thrusterImpulse(k), sprintf('%.2f', res.thrusterImpulse(k)), ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 7);
    end
end
title(ax, 'Per-thruster impulse (colour = pod)', 'FontWeight', 'normal');
end
