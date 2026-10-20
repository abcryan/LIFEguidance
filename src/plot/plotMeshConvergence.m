function fig = plotMeshConvergence(levels, p)
%PLOTMESHCONVERGENCE Total impulse and terminal errors per refinement level.
N = [levels.N];
J = arrayfun(@(l) l.sol.J, levels);
term = [levels.term];
fig = newFigure('mesh', [50 50 1300 550]);
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact');
ax = nexttile(tl); hold(ax, 'on'); grid(ax, 'on');
acc = [levels.accepted];
plot(ax, N, J, 'k-', 'LineWidth', 1.3);
plot(ax, N(acc), J(acc), 'ko', 'MarkerFaceColor', 'k');
plot(ax, N(~acc), J(~acc), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
for k = 1:numel(N)
    text(ax, N(k), J(k), sprintf('  L%d', levels(k).level), 'VerticalAlignment', 'bottom');
end
xlabel(ax, 'Number of intervals N'); ylabel(ax, 'Total impulse J (N s)');
title(ax, 'Total impulse vs. mesh (filled: accepted, red x: rejected)', 'FontWeight', 'normal');
ax = nexttile(tl); hold(ax, 'on'); grid(ax, 'on'); set(ax, 'YScale', 'log');
semilogy(ax, N, [term.r]/p.gate.r, 'o-', 'LineWidth', 1.3);
semilogy(ax, N, [term.v]/p.gate.v, 's-', 'LineWidth', 1.3);
semilogy(ax, N, [term.thetaAtan]/p.gate.theta, 'd-', 'LineWidth', 1.3);   % acos form floors at 0
semilogy(ax, N, [term.w]/p.gate.w, '^-', 'LineWidth', 1.3);
yline(ax, 1, 'k--', 'tolerance');
legend(ax, {'|r|/10^{-2} m', '|v|/10^{-5} m/s', '\theta/10^{-4} rad (atan2 form)', ...
            '|\omega|/10^{-7} s^{-1}'}, 'Location', 'southwest');
xlabel(ax, 'Number of intervals N'); ylabel(ax, 'Terminal error / tolerance');
title(ax, 'Propagated terminal errors (59)-(60)', 'FontWeight', 'normal');
end
