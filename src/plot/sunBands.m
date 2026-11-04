function sunBands(ax, p, tLim, intervals, s)
%SUNBANDS Field of regard of eq. (48) as a shaded band of anti-Sun angles
%   over tLim, its bounds, and red shading of the intervals outside it.
lo = 180 - rad2deg(p.sun.thetaMax); hi = 180 - rad2deg(p.sun.thetaMin);
for b = 1:numel(intervals)
    I = intervals(b);
    patch(ax, [I.t0 I.t1 I.t1 I.t0], [0 0 180 180], s.viol, 'EdgeColor', 'none', ...
          'HandleVisibility', 'off');
end
patch(ax, [tLim(1) tLim(2) tLim(2) tLim(1)], [lo lo hi hi], s.band, 'EdgeColor', 'none', ...
      'HandleVisibility', 'off');
yline(ax, [lo hi], 'g-', {sprintf('%.0f\\circ', lo), sprintf('%.0f\\circ', hi)}, 'LineWidth', 1.2, ...
      'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');
ylim(ax, [0 180]); yticks(ax, 0:30:180);
end
