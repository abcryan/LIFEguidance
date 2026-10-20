function videoFormation(res, p, sched, file)
%VIDEOFORMATION Animation in the virtual formation frame O^f: stationary
%   target point, approach path, spacecraft attitude, formation and body
%   axes, exhaust plumes.
mag = 3; plumeLen = 4;
x = res.prop.x;
rf = x(8:10,:) + p.rfp;                  % r^f_fi

fig = newFigure('videoFormation', [0 0 1280 720]);
pan = videoPanels(fig, res, p, 'Virtual formation frame O^f (rotating with the formation)');
ax = axes(fig, 'Position', [0.03 0.06 0.62 0.85]); hold(ax, 'on');
plot3(ax, rf(1,:), rf(2,:), rf(3,:), ':', 'Color', [0.4 0.4 0.4]);
trail = plot3(ax, nan, nan, nan, 'k', 'LineWidth', 1.5);
plot3(ax, p.rfp(1), p.rfp(2), p.rfp(3), 'o', 'MarkerSize', 10, ...
      'MarkerFaceColor', [0.2 0.4 0.9], 'Color', 'k');
text(ax, p.rfp(1), p.rfp(2), p.rfp(3) - 6, 'target point p', 'FontSize', 9);
lim = [min([rf, zeros(3,1)], [], 2) - 15, max([rf, zeros(3,1)], [], 2) + 15];
L = 0.25*(lim(1,2) - lim(1,1));
cols = [0.85 0.1 0.1; 0.1 0.6 0.1; 0.1 0.2 0.9];
nm = {'x_f', 'y_f', 'z_f'};
for i = 1:3
    e = L*((1:3)' == i);
    plot3(ax, [0 e(1)], [0 e(2)], [0 e(3)], '-', 'Color', cols(i,:), 'LineWidth', 2);
    text(ax, e(1), e(2), e(3), nm{i}, 'Color', cols(i,:), 'FontWeight', 'bold');
end
h = drawSpacecraft(ax, p, mag);
axis(ax, 'equal'); axis(ax, reshape(lim', 1, [])); grid(ax, 'on'); view(ax, 40, 22);
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
    [xs, T] = stateAt(res, t);
    pos = xs(8:10) + p.rfp;
    setSpacecraftPose(h, quatToDcm(xs(1:4)), pos, T, plumeLen);
    upto = res.prop.t <= t;
    set(trail, 'XData', [rf(1,upto) pos(1)], 'YData', [rf(2,upto) pos(2)], 'ZData', [rf(3,upto) pos(3)]);
    updatePanels(pan, t, T, sched.label(f), xs, p);
    writeVideo(vw, grabFrame(fig));
end
close(fig);
end
