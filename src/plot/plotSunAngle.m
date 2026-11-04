function fig = plotSunAngle(res, p)
%PLOTSUNANGLE Anti-Sun angle 180 - theta_i of n^i (propagated, exact Sun
%   direction) against the field of regard of eq. (48), with the intervals
%   outside it (red) and, for Problem 2, the points where it is active.
s = plotStyle(); sn = res.sun;
fig = newFigure('sunAngle', [50 50 1300 500]);
ax = axes(fig); hold(ax, 'on'); grid(ax, 'on');
sunBands(ax, p, [0 res.sol.tau], sn.intervals, s);
plot(ax, sn.t, sn.antiDeg, 'k', 'LineWidth', 1.5);
if ~isempty(sn.active)
    plot(ax, [sn.active.t], [sn.active.antiDeg], 'o', 'MarkerSize', 9, ...
         'MarkerFaceColor', s.prob(2,:), 'Color', 'k');
end
xlim(ax, [0 res.sol.tau]); xlabel(ax, 'Time t (s)');
ylabel(ax, 'Anti-Sun angle of z^i, 180\circ - \theta_i (deg)');
str = sprintf('max violation %.3f deg, %.0f s outside', sn.maxViolDeg, sn.timeOutside);
if res.problem >= 2
    str = sprintf('%s; active at %d points (dots)', str, numel(sn.active));
else
    str = [str, ' (not imposed in Problem 1)'];
end
title(ax, sprintf('Problem %d: Sun angle, eq. (48) — %s', res.problem, str), 'FontWeight', 'normal');
end
