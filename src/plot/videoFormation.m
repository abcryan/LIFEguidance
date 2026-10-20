function videoFormation(res, p, sched, file)
%VIDEOFORMATION Animation in the virtual formation frame O^f: stationary
%   virtual formation and target point, the formation's orbital
%   trajectory (a line through the origin), approach path, spacecraft
%   attitude, formation and body axes, exhaust plumes.
mag = 3; plumeLen = 4;
x = res.prop.x;
rf = x(8:10,:) + p.rfp;                  % r^f_fi

fig = newFigure('videoFormation', [50 50 1280 720]); fig.Resize = 'off';
pan = videoPanels(fig, res, p, 'Virtual formation frame O^f (rotating with the formation)');
ax = axes(fig, 'Position', [0.17 0.09 0.49 0.76]); hold(ax, 'on');
plot3(ax, rf(1,:), rf(2,:), rf(3,:), ':', 'Color', [0.4 0.4 0.4]);
trail = plot3(ax, nan, nan, nan, 'k', 'LineWidth', 1.5);
plot3(ax, p.rfp(1), p.rfp(2), p.rfp(3), 'o', 'MarkerSize', 10, ...
      'MarkerFaceColor', [0.2 0.4 0.9], 'Color', 'k');
text(ax, p.rfp(1), p.rfp(2), p.rfp(3) - 6, 'target point p', 'FontSize', 9);
drawFormation(ax, p);
ho = drawOrbitLine(ax);
V = [p.viz.collectors, p.viz.combiner];
lim = [min([rf, V], [], 2) - 15, max([rf, V], [], 2) + 15];
L = 0.25*(lim(1,2) - lim(1,1));
cols = [0.85 0.1 0.1; 0.1 0.6 0.1; 0.1 0.2 0.9];
nm = {'x_f', 'y_f', 'z_f'};
for i = 1:3
    e = L*((1:3)' == i);
    plot3(ax, [0 e(1)], [0 e(2)], [0 e(3)], '-', 'Color', cols(i,:), 'LineWidth', 2);
    text(ax, e(1), e(2), e(3), nm{i}, 'Color', cols(i,:), 'FontWeight', 'bold');
end
h = drawSpacecraft(ax, p, mag);
limV = reshape(lim', 1, []);
axis(ax, 'equal'); axis(ax, limV); grid(ax, 'on'); view(ax, 40, 22);
camlight(ax, 'headlight'); lighting(ax, 'gouraud');
xlabel(ax, 'x_f (m)'); ylabel(ax, 'y_f (m)'); zlabel(ax, 'z_f (m)');
title(ax, {'formation axes at origin; body axes red/green/blue', ...
    sprintf('spacecraft drawn x%d, plume %g m at T_{max}', mag, plumeLen*mag)}, 'FontSize', 9);

vw = VideoWriter(file, 'MPEG-4'); vw.FrameRate = sched.fps; vw.Quality = 75;
open(vw); cleanup = onCleanup(@() close(vw));
warnState = warning('off', 'MATLAB:hg:DiceyTransformMatrix');
restore = onCleanup(@() warning(warnState));
for f = 1:numel(sched.t)
    t = sched.t(f);
    xs = stateAt(res, t);
    T = thrustInWindow(res, sched.t(max(f - 1, 1)), t);   % peak since last frame
    pos = xs(8:10) + p.rfp;
    setSpacecraftPose(h, quatToDcm(xs(1:4)), pos, T, plumeLen);
    fr = frames(t, p);
    vb = interp1(res.prop.t, res.ref(4:6,:)', t)';          % formation velocity in O^b
    setOrbitLine(ho, fr.Rbf'*vb, norm(vb), limV);
    upto = res.prop.t <= t;
    set(trail, 'XData', [rf(1,upto) pos(1)], 'YData', [rf(2,upto) pos(2)], 'ZData', [rf(3,upto) pos(3)]);
    updatePanels(pan, t, T, sched.label, xs, p);
    writeVideo(vw, grabFrame(fig));
end
end
