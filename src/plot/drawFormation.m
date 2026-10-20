function h = drawFormation(ax, p)
%DRAWFORMATION Virtual formation pyramid: collectors (base corners) and
%   combiner (apex) as dots, dotted edges, dashed z_f axis from the base
%   centre to the apex. Place it with setFormationPose.
red = [0.8 0.1 0.1];
h.V = [p.viz.collectors, p.viz.combiner];            % 3x5, apex last
h.base = plot3(ax, nan, nan, nan, ':', 'Color', red, 'LineWidth', 1.5);
h.edges = plot3(ax, nan, nan, nan, ':', 'Color', red, 'LineWidth', 1.5);
h.axis = plot3(ax, nan, nan, nan, '--', 'Color', [0.55 0.05 0.05], 'LineWidth', 1.2);
h.dots = plot3(ax, nan, nan, nan, 'o', 'MarkerSize', 6, 'MarkerFaceColor', red, ...
               'MarkerEdgeColor', red, 'LineStyle', 'none');
setFormationPose(h, eye(3));
end
