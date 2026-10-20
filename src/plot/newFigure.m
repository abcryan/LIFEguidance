function fig = newFigure(name, pos)
%NEWFIGURE White figure; visible in the MATLAB desktop, invisible in -batch.
if nargin < 2, s = plotStyle(); pos = s.figPos; end
vis = 'off';
if feature('ShowFigureWindows'), vis = 'on'; end
fig = figure('Name', name, 'Color', 'w', 'Visible', vis, 'Position', pos);
end
