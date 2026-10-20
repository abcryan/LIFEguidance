function fig = newFigure(name, pos)
%NEWFIGURE Invisible white figure (works in -batch mode).
if nargin < 2, s = plotStyle(); pos = s.figPos; end
fig = figure('Name', name, 'Color', 'w', 'Visible', 'off', 'Position', pos);
end
