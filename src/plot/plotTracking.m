function fig = plotTracking(res, p)
%PLOTTRACKING Eigenaxis attitude error and components of w^f_fi, r^f_pi,
%   rdot^f_pi (propagated), with magnified arrival views of the norms
%   against the terminal tolerances (59)-(60).
s = plotStyle();
t = res.prop.t; x = res.prop.x; tau = res.sol.tau;
fig = newFigure('tracking');
tl = tiledlayout(fig, 4, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('Tracking errors (independent propagation), \\tau = %.0f s', tau));
data = {rad2deg(res.thetaAcos), '\theta_i (deg)', res.thetaAcos, p.gate.theta, 'rad'
        x(5:7,:)*1e3, '\omega^f_{fi} (mrad/s)', vecnorm(x(5:7,:)), p.gate.w, '1/s'
        x(8:10,:), 'r^f_{pi} (m)', vecnorm(x(8:10,:)), p.gate.r, 'm'
        x(11:13,:)*1e3, 'dr^f_{pi}/dt (mm/s)', vecnorm(x(11:13,:)), p.gate.v, 'm/s'};
zoom = t >= tau - 0.1*tau;
for i = 1:4
    ax = nexttile(tl, [1 2]); hold(ax, 'on'); grid(ax, 'on');
    Y = data{i, 1};
    if size(Y, 1) == 1
        plot(ax, t, Y, 'k', 'LineWidth', 1.2);
    else
        for c = 1:3
            plot(ax, t, Y(c,:), 'Color', s.comp(c,:), 'LineWidth', 1.2);
        end
        legend(ax, {'x', 'y', 'z'}, 'Location', 'northeast');
    end
    ylabel(ax, data{i, 2}); xlim(ax, [0 tau]);
    if i == 4, xlabel(ax, 'Time t (s)'); end

    ax = nexttile(tl); hold(ax, 'on'); grid(ax, 'on');
    nrm = max(data{i, 3}(zoom), realmin);
    patch(ax, [t(find(zoom, 1)) tau tau t(find(zoom, 1))], [1e-12 1e-12 data{i,4} data{i,4}], ...
          s.band, 'EdgeColor', 'none');
    semilogy(ax, t(zoom), nrm, 'k', 'LineWidth', 1.2);
    yline(ax, data{i, 4}, 'g-', 'tolerance', 'LineWidth', 1.2, 'LabelHorizontalAlignment', 'left');
    plot(ax, tau, max(data{i, 3}(end), realmin), 'ro', 'MarkerFaceColor', 'r');
    set(ax, 'YScale', 'log');
    ylim(ax, [min(data{i, 4}/100, max(min(nrm), 1e-12)), max(nrm)*2]);
    xlim(ax, [t(find(zoom, 1)) tau]);
    ylabel(ax, sprintf('norm (%s)', data{i, 5}));
    title(ax, sprintf('arrival: %.2e %s', data{i, 3}(end), data{i, 5}), 'FontWeight', 'normal');
    if i == 4, xlabel(ax, 'Time t (s)'); end
end
end
