function drawHalo(ax, res, p)
%DRAWHALO Reference L2 halo orbit in the SER frame (relative to Earth)
%   with the arc flown during the transfer highlighted.
hold(ax, 'on');
rel = @(r) (r - p.rbe/p.dse)*p.dse/1e6;
H = rel(res.haloFull(1:3,:)); A = rel(res.ref(1:3,:)/p.dse);
plot3(ax, H(1,:), H(2,:), H(3,:), 'Color', [0.6 0.6 0.6]);
plot3(ax, A(1,:), A(2,:), A(3,:), 'r', 'LineWidth', 3);
plot3(ax, 0, 0, 0, 'bo', 'MarkerFaceColor', 'b'); text(ax, 0, 0, 0, '  Earth');
axis(ax, 'equal'); grid(ax, 'on'); view(ax, 35, 25);
xlabel(ax, 'x_b (10^3 km)'); ylabel(ax, 'y_b (10^3 km)'); zlabel(ax, 'z_b (10^3 km)');
title(ax, sprintf('L_2 halo (period %.2f d), red: %.1f h arc flown', p.haloPeriodDays, res.sol.tau/3600));
end
