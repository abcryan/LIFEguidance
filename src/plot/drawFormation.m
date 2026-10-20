function h = drawFormation(ax, V, col)
%DRAWFORMATION Virtual formation pyramid with points V (3x5, corners then
%   apex, in the formation's own frame): collectors and combiner as dots,
%   dotted edges, dashed axis from the base centre to the apex, in colour
%   col (default red). Place it with setFormationPose.
if nargin < 3, col = [0.8 0.1 0.1]; end
h.V = V;
h.base = plot3(ax, nan, nan, nan, ':', 'Color', col, 'LineWidth', 1.5);
h.edges = plot3(ax, nan, nan, nan, ':', 'Color', col, 'LineWidth', 1.5);
h.axis = plot3(ax, nan, nan, nan, '--', 'Color', 0.7*col, 'LineWidth', 1.2);
h.dots = plot3(ax, nan, nan, nan, 'o', 'MarkerSize', 6, 'MarkerFaceColor', col, ...
               'MarkerEdgeColor', col, 'LineStyle', 'none');
setFormationPose(h, eye(3));
end
