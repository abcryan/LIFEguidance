function figs = plotComparison(R, p, outDir)
%PLOTCOMPARISON Figures comparing the solutions in the cell array R (one
%   validated result per problem), saved as PNG to outDir.
list = {
    'sun_angle_comparison',  @() sunAngleFig(R, p)
    'sky_plot',              @() skyFig(R, p)
    'tracking_comparison',   @() trackingFig(R, p)
    'impulse_comparison',    @() impulseFig(R, p)
    'trajectory_comparison', @() trajectoryFig(R, p)
    };
for k = 1:size(list, 1)
    figs.(list{k, 1}) = list{k, 2}();
    exportgraphics(figs.(list{k, 1}), fullfile(outDir, [list{k, 1} '.png']), 'Resolution', 200);
end
end

function fig = sunAngleFig(R, p)
s = plotStyle();
fig = newFigure('sunAngleCompare', [50 50 1300 550]);
ax = axes(fig); hold(ax, 'on'); grid(ax, 'on');
tMax = max(cellfun(@(r) r.sol.tau, R));
sunBands(ax, p, [0 tMax], R{1}.sun.intervals, s);   % red: Problem 1 outside
h = gobjects(1, numel(R)); lbl = cell(1, numel(R));
for i = 1:numel(R)
    sn = R{i}.sun; c = s.prob(R{i}.problem,:);
    h(i) = plot(ax, sn.t, sn.antiDeg, '-', 'Color', c, 'LineWidth', 1.8);
    if ~isempty(sn.active)
        plot(ax, [sn.active.t], [sn.active.antiDeg], 'o', 'MarkerSize', 9, 'MarkerFaceColor', c, ...
             'Color', 'k', 'HandleVisibility', 'off');
    end
    lbl{i} = sprintf('Problem %d: max violation %.3f\\circ, %.0f s outside, active at %d points', ...
                     R{i}.problem, sn.maxViolDeg, sn.timeOutside, numel(sn.active));
end
legend(ax, h, lbl, 'Location', 'northwest');
xlim(ax, [0 tMax]); xlabel(ax, 'Time t (s)');
ylabel(ax, 'Anti-Sun angle of z^i, 180\circ - \theta_i (deg)');
title(ax, {'Sun angle along the propagated trajectories (green: field of regard, eq. (45))', ...
    sprintf('red: Problem %d outside the field of regard; dots: constraint active', R{1}.problem)}, ...
    'FontWeight', 'normal');
end

function fig = skyFig(R, p)
s = plotStyle();
fig = newFigure('sky', [50 50 800 800]);
ax = axes(fig); drawSky(ax, p);
h = gobjects(1, numel(R));
for i = 1:numel(R)
    P = skyXY(R{i}.sun); c = s.prob(R{i}.problem,:);
    h(i) = plot(ax, P(1,:), P(2,:), '-', 'Color', c, 'LineWidth', 1.8);
    plot(ax, P(1,1), P(2,1), 's', 'MarkerFaceColor', c, 'Color', 'k', 'MarkerSize', 8);
    plot(ax, P(1,end), P(2,end), 'o', 'MarkerFaceColor', c, 'Color', 'k', 'MarkerSize', 8);
    for tk = 0:1800:R{i}.sol.tau             % half-hour ticks
        [~, k] = min(abs(R{i}.sun.t - tk));
        plot(ax, P(1,k), P(2,k), '.', 'Color', c, 'MarkerSize', 14, 'HandleVisibility', 'off');
    end
end
legend(ax, h, arrayfun(@(i) sprintf('Problem %d', R{i}.problem), 1:numel(R), 'UniformOutput', false), ...
       'Location', 'southwest');
title(ax, {'Direction of z^i about the Sun line (radius: anti-Sun angle)', ...
    'square: start, circle: arrival, dots every 30 min; green: field of regard'}, 'FontWeight', 'normal');
end

function fig = trackingFig(R, p)
s = plotStyle();
fig = newFigure('trackingCompare', [50 50 1300 800]);
tl = tiledlayout(fig, 3, 1, 'TileSpacing', 'compact');
yl = {'Attitude error \xi_i (deg)', '|r^f_{pi}| (m)', '|\omega^f_{fi}| (mrad/s)'};
for j = 1:3
    ax = nexttile(tl); hold(ax, 'on'); grid(ax, 'on');
    for i = 1:numel(R)
        r = R{i}; x = r.prop.x;
        y = {rad2deg(r.thetaAtan), vecnorm(x(8:10,:)), 1e3*vecnorm(x(5:7,:))};
        plot(ax, r.prop.t, y{j}, 'Color', s.prob(r.problem,:), 'LineWidth', 1.5, ...
             'DisplayName', sprintf('Problem %d', r.problem));
    end
    ylabel(ax, yl{j});
    if j == 1, legend(ax, 'Location', 'northeast'); end
end
xlabel(ax, 'Time t (s)');
title(tl, sprintf('Tracking errors (independent propagation); arrival gates %g rad, %g m, %g 1/s', ...
      p.gate.theta, p.gate.r, p.gate.w));
end

function fig = impulseFig(R, p)
s = plotStyle();
fig = newFigure('impulseCompare', [50 50 1300 550]);
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact');
ax = nexttile(tl); hold(ax, 'on'); grid(ax, 'on');
for i = 1:numel(R)
    r = R{i};
    plot(ax, r.sol.tNode, r.cumImpulse, 'Color', s.prob(r.problem,:), 'LineWidth', 1.8, ...
         'DisplayName', sprintf('Problem %d: %.2f N s (\\Deltam %.3f kg)', r.problem, r.J, r.dm));
end
legend(ax, 'Location', 'northwest'); xlabel(ax, 'Time t (s)'); ylabel(ax, 'Cumulative impulse (N s)');
ax = nexttile(tl); hold(ax, 'on'); grid(ax, 'on');
I = cell2mat(cellfun(@(r) r.thrusterImpulse(:), R, 'UniformOutput', false));
bh = bar(ax, 1:p.nThr, I, 'grouped');
for i = 1:numel(R)
    bh(i).FaceColor = s.prob(R{i}.problem,:);
    bh(i).DisplayName = sprintf('Problem %d', R{i}.problem);
end
set(ax, 'XTick', 1:p.nThr, 'XTickLabel', p.thr.label, 'XTickLabelRotation', 60);
ylabel(ax, 'Integrated impulse (N s)'); legend(ax, 'Location', 'northeast');
title(tl, 'Impulse');
end

function fig = trajectoryFig(R, p)
% In the prescribed formation frame (that of p): a solution with formation
% phase phase0_i is rotated by Rz(phase0_i - p.target.phase0) into it, so
% the target points of different arrival phases appear on their circle
% about z_f (Sec. 5.3, Problem 3)
s = plotStyle();
fig = newFigure('trajCompare', [50 50 1100 850]);
ax = axes(fig); hold(ax, 'on'); grid(ax, 'on');
nR = numel(R); h = gobjects(1, nR); hp = gobjects(1, nR); lblP = cell(1, nR);
th = linspace(0, 2*pi, 181); rr = norm(p.rfp(1:2));
plot3(ax, rr*cos(th), rr*sin(th), p.rfp(3)*ones(size(th)), ':', 'Color', [0.5 0.5 0.5], ...
      'HandleVisibility', 'off');
for i = 1:nR
    r = R{i}; c = s.prob(r.problem,:);
    Rz = rot3(r.p.target.phase0 - p.target.phase0);
    rf = Rz*(r.prop.x(8:10,:) + p.rfp);
    h(i) = plot3(ax, rf(1,:), rf(2,:), rf(3,:), 'Color', c, 'LineWidth', 1.8);
    k = round(linspace(1, numel(r.prop.t), 13));
    Z = zeros(3, numel(k));
    for j = 1:numel(k)
        Z(:,j) = Rz*quatToDcm(r.prop.x(1:4,k(j)))*p.sun.n;
    end
    Z = 0.12*rr*Z;                              % arrow length [m]
    quiver3(ax, rf(1,k), rf(2,k), rf(3,k), Z(1,:), Z(2,:), Z(3,:), 0, ...
            'Color', c, 'LineWidth', 1.2, 'MaxHeadSize', 0.4, 'HandleVisibility', 'off');
    tp = Rz*p.rfp;
    hp(i) = plot3(ax, tp(1), tp(2), tp(3), 'o', 'MarkerSize', 11, 'MarkerFaceColor', c, 'Color', 'k');
    lblP{i} = sprintf('target point, Problem %d (\\phi_\\tau = %.1f\\circ)', r.problem, ...
                      rad2deg(angle(exp(1i*r.sol.phiTau))));
end
hs = plot3(ax, rf(1,1), rf(2,1), rf(3,1), 'ks', 'MarkerFaceColor', 'g', 'MarkerSize', 9);
cols = [0.85 0.1 0.1; 0.1 0.6 0.1; 0.1 0.2 0.9];
for j = 1:3
    e = 0.3*rr*((1:3)' == j);
    plot3(ax, [0 e(1)], [0 e(2)], [0 e(3)], 'Color', cols(j,:), 'LineWidth', 2, 'HandleVisibility', 'off');
end
axis(ax, 'equal'); view(ax, 40, 22);
xlabel(ax, 'x_f (m)'); ylabel(ax, 'y_f (m)'); zlabel(ax, 'z_f (m)');
legend(ax, [h, hs, hp], ...
       [arrayfun(@(i) sprintf('Problem %d', R{i}.problem), 1:nR, 'UniformOutput', false), ...
        {'start'}, lblP], 'Location', 'southoutside', 'NumColumns', 2);
title(ax, {'Approach in the prescribed formation frame O^f (\phi(0) = \phi_0); dotted: circle of the target point', ...
    'arrows: body axis z^i every \tau/12; formation axes at the origin'}, 'FontWeight', 'normal');
end
