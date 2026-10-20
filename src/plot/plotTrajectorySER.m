function fig = plotTrajectorySER(res, p)
%PLOTTRAJECTORYSER Halo reference and the formation-centred close-up in
%   SER-parallel axes, with attitude snapshots.
mag = 3;
x = res.prop.x; nP = numel(res.prop.t);
dB = zeros(3, nP); tgtB = zeros(3, nP);
for k = 1:nP
    dB(:,k) = res.Rbf(:,:,k)*(x(8:10,k) + p.rfp);
    tgtB(:,k) = res.Rbf(:,:,k)*p.rfp;
end
fig = newFigure('trajSER', [50 50 1500 750]);
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact');
ax = nexttile(tl); hold(ax, 'on');
rel = @(r) (r - p.rbe/p.dse)*p.dse/1e6;
H = rel(res.haloFull(1:3,:)); A = rel(res.ref(1:3,:)/p.dse);
plot3(ax, H(1,:), H(2,:), H(3,:), 'Color', [0.6 0.6 0.6]);
plot3(ax, A(1,:), A(2,:), A(3,:), 'r', 'LineWidth', 3);
plot3(ax, 0, 0, 0, 'bo', 'MarkerFaceColor', 'b'); text(ax, 0, 0, 0, '  Earth');
axis(ax, 'equal'); grid(ax, 'on'); view(ax, 35, 25);
xlabel(ax, 'x_b (10^3 km)'); ylabel(ax, 'y_b (10^3 km)'); zlabel(ax, 'z_b (10^3 km)');
title(ax, sprintf('L_2 halo (period %.2f d), red: %.1f h arc flown', p.haloPeriodDays, res.sol.tau/3600));

ax = nexttile(tl); hold(ax, 'on');
plot3(ax, dB(1,:), dB(2,:), dB(3,:), 'k', 'LineWidth', 1.5);
plot3(ax, tgtB(1,:), tgtB(2,:), tgtB(3,:), 'Color', [0.2 0.4 0.9], 'LineWidth', 2);
plot3(ax, 0, 0, 0, 'k+', 'MarkerSize', 10, 'LineWidth', 2);
for t = linspace(0, res.sol.tau, 7)
    [xs, T] = stateAt(res, t);
    fr = frames(t, p);
    h = drawSpacecraft(ax, p, mag);
    setSpacecraftPose(h, fr.Rbf*quatToDcm(xs(1:4)), fr.Rbf*(xs(8:10) + p.rfp), T, 4);
    tp = fr.Rbf*p.rfp;
    plot3(ax, tp(1), tp(2), tp(3), 'o', 'MarkerFaceColor', [0.2 0.4 0.9], 'Color', 'k');
end
axis(ax, 'equal'); grid(ax, 'on'); view(ax, 35, 20);
camlight(ax, 'headlight'); lighting(ax, 'gouraud');
xlabel(ax, 'x_b (m)'); ylabel(ax, 'y_b (m)'); zlabel(ax, 'z_b (m)');
legend(ax, {'spacecraft path', 'target point path', 'formation centre'}, 'Location', 'northeast');
title(ax, sprintf('Formation-centred close-up, SER-parallel axes (spacecraft x%d)', mag));
end
