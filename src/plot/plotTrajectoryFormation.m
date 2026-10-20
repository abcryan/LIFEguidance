function fig = plotTrajectoryFormation(res, p)
%PLOTTRAJECTORYFORMATION Approach path in the formation frame with
%   attitude snapshots (spacecraft drawn magnified).
mag = 3;
x = res.prop.x; rf = x(8:10,:) + p.rfp;
fig = newFigure('trajF', [50 50 1100 850]);
ax = axes(fig); hold(ax, 'on');
plot3(ax, rf(1,:), rf(2,:), rf(3,:), 'k', 'LineWidth', 1.5);
plot3(ax, p.rfp(1), p.rfp(2), p.rfp(3), 'o', 'MarkerSize', 10, 'MarkerFaceColor', [0.2 0.4 0.9], 'Color', 'k');
plot3(ax, rf(1,1), rf(2,1), rf(3,1), 'ks', 'MarkerFaceColor', 'g', 'MarkerSize', 9);
cols = [0.85 0.1 0.1; 0.1 0.6 0.1; 0.1 0.2 0.9];
for i = 1:3
    e = 25*((1:3)' == i);
    plot3(ax, [0 e(1)], [0 e(2)], [0 e(3)], 'Color', cols(i,:), 'LineWidth', 2);
end
snapT = linspace(0, res.sol.tau, 7);
for t = snapT
    [xs, T] = stateAt(res, t);
    h = drawSpacecraft(ax, p, mag);
    setSpacecraftPose(h, quatToDcm(xs(1:4)), xs(8:10) + p.rfp, T, 4);
    text(ax, xs(8) + p.rfp(1), xs(9), xs(10) + 14, sprintf('%.0f s', t), 'FontSize', 8);
end
axis(ax, 'equal'); grid(ax, 'on'); view(ax, 40, 22);
camlight(ax, 'headlight'); lighting(ax, 'gouraud');
xlabel(ax, 'x_f (m)'); ylabel(ax, 'y_f (m)'); zlabel(ax, 'z_f (m)');
legend(ax, {'spacecraft path', 'target point p', 'start'}, 'Location', 'northeast');
title(ax, sprintf('Approach in the formation frame O^f (spacecraft drawn x%d, axes: x_f red, y_f green, z_f blue)', mag));
end
