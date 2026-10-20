function h = drawOrbitLine(ax)
%DRAWORBITLINE Orbital (halo) trajectory of the formation centre as seen in
%   a formation-centred view: at true scale over the plotted region it is a
%   straight line through the origin along the direction of motion.
%   Place it with setOrbitLine.
col = [0.85 0.33 0.1];
h.line = plot3(ax, nan, nan, nan, '-', 'Color', col, 'LineWidth', 1.8);
h.label = text(ax, nan, nan, nan, '', 'Color', col, 'FontSize', 9, 'FontWeight', 'bold', ...
               'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
end
