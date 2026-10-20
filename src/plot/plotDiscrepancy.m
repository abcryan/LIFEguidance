function fig = plotDiscrepancy(res, p)
%PLOTDISCREPANCY Independent propagation minus collocation at the mesh
%   nodes, and the absolute-model cross-check.
tN = res.prop.t(res.prop.nodeIdx);
fig = newFigure('discrepancy', [50 50 1300 800]);
tl = tiledlayout(fig, 4, 1, 'TileSpacing', 'compact');
title(tl, 'Discrepancy of the independent propagation');
f = {'r', '|\Deltar| (m)', p.gate.r; 'v', '|\Deltav| (m/s)', p.gate.v;
     'theta', '\Delta\theta (rad)', p.gate.theta; 'w', '|\Delta\omega| (1/s)', p.gate.w};
for i = 1:4
    ax = nexttile(tl); hold(ax, 'on'); grid(ax, 'on'); set(ax, 'YScale', 'log');
    semilogy(ax, tN, max(res.disc.(f{i,1}), 1e-20), 'b', 'LineWidth', 1.2);
    semilogy(ax, tN, max(res.absDisc.(f{i,1}), 1e-20), 'Color', [0.85 0.33 0.1], 'LineWidth', 1.2);
    yline(ax, f{i,3}, 'g-', 'arrival tolerance');
    yline(ax, p.gate.absRel*f{i,3}, 'g--', 'absolute-model gate');
    ylabel(ax, f{i,2}); xlim(ax, [0 res.sol.tau]);
    if i == 1
        legend(ax, {'propagation - collocation', 'absolute model - relative propagation'}, ...
               'Location', 'southeast');
    end
end
xlabel(ax, 'Time t (s)');
end
