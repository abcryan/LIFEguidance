function fig = plotHaloOrbit(res, p)
%PLOTHALOORBIT Where the virtual formation flies: the reference L2 halo
%   orbit in the Sun-Earth rotating frame O^b (centred on L2), with Earth,
%   the Sun direction, the arc of the transfer, the formation axes
%   x_f, y_f, z_f and the line of sight z_f to the target star; two
%   projections; and a zoom on the arc showing the formation spin.
axCol = [0.85 0.1 0.1; 0.1 0.6 0.1; 0.1 0.2 0.9];  % x_f, y_f, z_f as in the videos
xL2 = l2Location(p.eta);
km = p.dse/1e6;                                   % normalized -> 10^3 km
H = (res.haloFull(1:3,:) - [xL2; 0; 0])*km;
A = (res.ref(1:3,:)/p.dse - [xL2; 0; 0])*km;
E = ([1 - p.eta; 0; 0] - [xL2; 0; 0])*km;         % Earth
period = p.haloPeriodDays;

fig = newFigure('halo', [50 50 1400 950]);
tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, {sprintf(['Reference L_2 halo orbit of the virtual formation (period %.2f d, ' ...
    'SER frame O^b centred on L_2)'], period), ...
    sprintf(['yellow: formation at the start of the transfer; the %.1f h arc flown (black) ' ...
    'is too short to see at orbit scale, see the zoom'], res.sol.tau/3600)});

% 3D overview with the formation frame at the start of the transfer
ax = nexttile(tl); hold(ax, 'on'); grid(ax, 'on');
haloAndArc(ax, H, A, [1 2 3]);
plot3(ax, 0, 0, 0, 'k+', 'MarkerSize', 10, 'LineWidth', 2);
text(ax, 0, 0, 0, '  L_2', 'FontWeight', 'bold');
Lax = 0.35*max(range(H, 2));
fr = frames(0, p);
nm = {'x_f', 'y_f', 'z_f'};
for i = 1:3
    e = Lax*fr.Rbf(:,i);
    quiver3(ax, A(1,1), A(2,1), A(3,1), e(1), e(2), e(3), 0, 'Color', axCol(i,:), ...
            'LineWidth', 2, 'MaxHeadSize', 0.4);
    text(ax, A(1,1) + 1.08*e(1), A(2,1) + 1.08*e(2), A(3,1) + 1.08*e(3), nm{i}, ...
         'Color', axCol(i,:), 'FontWeight', 'bold');
end
los = 2.2*Lax*fr.Rbf(:,3);
plot3(ax, A(1,1) + [0 los(1)], A(2,1) + [0 los(2)], A(3,1) + [0 los(3)], ':', ...
      'Color', axCol(3,:), 'LineWidth', 1.2);
text(ax, A(1,1) + los(1), A(2,1) + los(2), A(3,1) + los(3), '  to target star (LOS)', ...
     'Color', axCol(3,:));
sunArrow(ax, H, true);
axis(ax, 'equal'); view(ax, 40, 22);
xlabel(ax, 'x_b (10^3 km)'); ylabel(ax, 'y_b (10^3 km)'); zlabel(ax, 'z_b (10^3 km)');
title(ax, 'Halo orbit and formation axes at t = 0 (axes not to scale)', 'FontWeight', 'normal');

% Projections, with Earth in the ecliptic view
ax = nexttile(tl); hold(ax, 'on'); grid(ax, 'on');
haloAndArc(ax, H, A, [1 2]);
plot(ax, 0, 0, 'k+', 'MarkerSize', 10, 'LineWidth', 2); text(ax, 0, 0, '  L_2');
plot(ax, E(1), E(2), 'o', 'MarkerFaceColor', [0.1 0.3 0.9], 'Color', 'k', 'MarkerSize', 8);
text(ax, E(1), E(2), '  Earth');
sunArrow(ax, H, false);
axis(ax, 'equal'); xlabel(ax, 'x_b (10^3 km)'); ylabel(ax, 'y_b (10^3 km)');
title(ax, 'Ecliptic view (x_b-y_b), Earth and Sun direction', 'FontWeight', 'normal');

ax = nexttile(tl); hold(ax, 'on'); grid(ax, 'on');
haloAndArc(ax, H, A, [2 3]);
plot(ax, 0, 0, 'k+', 'MarkerSize', 10, 'LineWidth', 2); text(ax, 0, 0, '  L_2');
axis(ax, 'equal'); xlabel(ax, 'y_b (10^3 km)'); ylabel(ax, 'z_b (10^3 km)');
title(ax, 'Seen from Earth (y_b-z_b)', 'FontWeight', 'normal');

% Zoom on the arc: formation axes at several times (km, relative to start)
ax = nexttile(tl); hold(ax, 'on'); grid(ax, 'on');
Z = (res.ref(1:3,:) - res.ref(1:3,1))/1e3;       % km
plot3(ax, Z(1,:), Z(2,:), Z(3,:), 'k', 'LineWidth', 2);
Lz = 0.18*max(vecnorm(Z));
for t = linspace(0, res.sol.tau, 5)
    k = find(res.prop.t >= t, 1); if isempty(k), k = numel(res.prop.t); end
    c = Z(:,k); R = res.Rbf(:,:,k);
    for i = 1:3
        e = Lz*R(:,i);
        plot3(ax, c(1) + [0 e(1)], c(2) + [0 e(2)], c(3) + [0 e(3)], 'Color', axCol(i,:), ...
              'LineWidth', 2);
    end
    plot3(ax, c(1), c(2), c(3), 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 5);
    text(ax, c(1), c(2), c(3) - 0.35*Lz, sprintf('%.0f h', res.prop.t(k)/3600), ...
         'HorizontalAlignment', 'center', 'FontSize', 9);
end
axis(ax, 'equal'); view(ax, 40, 22);
xlabel(ax, '\Deltax_b (km)'); ylabel(ax, '\Deltay_b (km)'); zlabel(ax, '\Deltaz_b (km)');
title(ax, sprintf(['Arc flown in %.1f h (%.0f km): formation axes x_f (red), y_f (green), ' ...
    'z_f = LOS (blue), spinning about z_f once per 24 h'], res.sol.tau/3600, ...
    sum(vecnorm(diff(res.ref(1:3,:), 1, 2)))/1e3), 'FontWeight', 'normal');
end

function haloAndArc(ax, H, A, dims)
if numel(dims) == 3
    plot3(ax, H(1,:), H(2,:), H(3,:), 'Color', [0.55 0.55 0.55], 'LineWidth', 1);
    plot3(ax, A(1,:), A(2,:), A(3,:), 'k', 'LineWidth', 4);
    plot3(ax, A(1,1), A(2,1), A(3,1), 'ko', 'MarkerFaceColor', 'y', 'MarkerSize', 8);
else
    plot(ax, H(dims(1),:), H(dims(2),:), 'Color', [0.55 0.55 0.55], 'LineWidth', 1);
    plot(ax, A(dims(1),:), A(dims(2),:), 'k', 'LineWidth', 4);
    plot(ax, A(dims(1),1), A(dims(2),1), 'ko', 'MarkerFaceColor', 'y', 'MarkerSize', 8);
end
end

function sunArrow(ax, H, is3D)
% The Sun lies along -x_b (far outside the plot)
L = 0.25*range(H(1,:)) + 0.15*max(range(H, 2));
x0 = min(H(1,:)) - 0.05*L;
if is3D
    quiver3(ax, x0, 0, 0, -L, 0, 0, 0, 'Color', [0.95 0.65 0], 'LineWidth', 2.5, 'MaxHeadSize', 0.6);
    text(ax, x0 - L, 0, 0, 'to Sun  ', 'HorizontalAlignment', 'right', 'Color', [0.85 0.5 0]);
else
    quiver(ax, x0, 0, -L, 0, 0, 'Color', [0.95 0.65 0], 'LineWidth', 2.5, 'MaxHeadSize', 0.6);
    text(ax, x0 - L, 0, 'to Sun  ', 'HorizontalAlignment', 'right', 'Color', [0.85 0.5 0]);
end
end

function x = l2Location(eta)
% Collinear point beyond Earth: zero of the x-gradient of the potential
f = @(x) x - (1 - eta)*(x + eta)/abs(x + eta)^3 - eta*(x - 1 + eta)/abs(x - 1 + eta)^3;
x = fzero(f, [1 - eta + 1e-4, 1.1]);
end
