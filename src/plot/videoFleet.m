function videoFleet(F, p, sched, file)
%VIDEOFLEET Problem 4 (or 5, F.problem) animation of all five spacecraft (Sec. 5.3), one
%   colour per spacecraft (plotStyle s.sc). Left: formation-centred view
%   with SER-parallel axes: the target formation spinning with the shared
%   phase (red) and its five moving target points, the initial formation
%   at t = 0 (grey, fixed), every path, and the spacecraft drawn magnified
%   with their plumes. Right, top: the view fixed to the target formation
%   frame (target points fixed, approach paths, z^i arrows). Right, middle:
%   thruster activity of all 5 x 20 thrusters (peak T/Tmax since the
%   previous frame). Right, bottom: distance to the target point, anti-Sun
%   angle of z^i and cumulative impulse, with a time cursor.
s = plotStyle(); nS = numel(F.res); R = F.res;
mag = 12; plumeLen = 4;
tau = R{1}.sol.tau;
uMax = max(cellfun(@(r) max(r.sol.U(:)), R));
lo = 180 - rad2deg(p.sun.thetaMax); hi = 180 - rad2deg(p.sun.thetaMin);

% Paths in SER-parallel axes and in the target formation frame
dB = cell(1, nS); tgtB = cell(1, nS); rF = cell(1, nS);
for i = 1:nS
    r = R{i}; n = numel(r.prop.t);
    dB{i} = zeros(3, n); tgtB{i} = zeros(3, n);
    for k = 1:n
        dB{i}(:,k) = r.Rbf(:,:,k)*(r.prop.x(8:10,k) + r.p.rfp);
        tgtB{i}(:,k) = r.Rbf(:,:,k)*r.p.rfp;
    end
    rF{i} = r.prop.x(8:10,:) + r.p.rfp;
end

fig = newFigure('videoFleet', [50 50 1280 720]); fig.Resize = 'off';
pn = 4; if isfield(F, 'problem'), pn = F.problem; end
extra = ''; if pn == 5, extra = sprintf(', collision avoidance d_{min} = %.0f m', p.collision.dmin); end
annotation(fig, 'textbox', [0.01 0.935 0.98 0.06], 'String', sprintf(['Problem %d: five spacecraft, ' ...
    'shared arrival time \\tau^* = %.0f s and phase \\phi_\\tau^* = %.2f\\circ, J_%d^* = %.2f N s%s'], ...
    pn, tau, rad2deg(R{1}.sol.phiTau), pn, sum(cellfun(@(r) r.J, R)), extra), ...
    'EdgeColor', 'none', 'FontSize', 15, 'FontWeight', 'bold');
hTime = annotation(fig, 'textbox', [0.01 0.885 0.60 0.05], 'EdgeColor', 'none', ...
    'FontSize', 12, 'FontName', 'FixedWidth');

% Main view: SER-parallel, formation-centred
ax = axes(fig, 'Position', [0.03 0.09 0.55 0.78]); hold(ax, 'on');
pk = R{1}.p;
hf = drawFormation(ax, [pk.viz.collectors, pk.viz.combiner]);
h0 = drawFormation(ax, initialFormationPoints(pk), [0.55 0.55 0.55]);
setFormationPose(h0, R{1}.Rbf(:,:,1));
plot3(ax, 0, 0, 0, 'k+', 'MarkerSize', 10, 'LineWidth', 2);
trail = gobjects(1, nS); tgt = gobjects(1, nS); sc = cell(1, nS); lbl = gobjects(1, nS);
for i = 1:nS
    c = s.sc(i,:);
    plot3(ax, dB{i}(1,:), dB{i}(2,:), dB{i}(3,:), ':', 'Color', 0.4 + 0.6*c);
    trail(i) = plot3(ax, nan, nan, nan, '-', 'Color', c, 'LineWidth', 2);
    tgt(i) = plot3(ax, nan, nan, nan, 'p', 'MarkerSize', 12, 'MarkerFaceColor', 'w', ...
                   'Color', c, 'LineWidth', 1.5);
    sc{i} = drawSpacecraft(ax, p, mag);
    lbl(i) = text(ax, 0, 0, 0, sprintf('%d', i), 'Color', c, 'FontWeight', 'bold', 'FontSize', 11, ...
                  'BackgroundColor', 'w', 'EdgeColor', c, 'Margin', 1, 'HorizontalAlignment', 'center');
end
lim = zoomLimits([cell2mat(dB), cell2mat(tgtB)], 60, 0.3);
axis(ax, 'equal'); axis(ax, lim); grid(ax, 'on'); view(ax, 35, 20);
camlight(ax, 'headlight'); lighting(ax, 'gouraud');
xlabel(ax, 'x_b (m)'); ylabel(ax, 'y_b (m)'); zlabel(ax, 'z_b (m)'); set(ax, 'FontSize', 8);
title(ax, sprintf(['SER-parallel axes, origin at the formation centre; spacecraft drawn x%d, ' ...
    'plume %g m at T_{max}; stars: target points'], mag, plumeLen*mag), 'FontSize', 9, 'FontWeight', 'normal');
legend(ax, [hf.dots, h0.dots, trail], [{'target formation', 'initial formation (t = 0)'}, ...
    arrayfun(@(i) sprintf('%d %s', i, p.fleet.name(i)), 1:nS, 'UniformOutput', false)], ...
    'Location', 'southoutside', 'NumColumns', 4, 'FontSize', 7, 'AutoUpdate', 'off');

% Formation-frame view
axF = axes(fig, 'Position', [0.62 0.56 0.36 0.33]); hold(axF, 'on');
drawFormation(axF, [pk.viz.collectors, pk.viz.combiner]);
dotF = gobjects(1, nS); zF = gobjects(1, nS); trF = gobjects(1, nS);
for i = 1:nS
    c = s.sc(i,:);
    plot3(axF, rF{i}(1,:), rF{i}(2,:), rF{i}(3,:), ':', 'Color', 0.4 + 0.6*c);
    trF(i) = plot3(axF, nan, nan, nan, '-', 'Color', c, 'LineWidth', 1.5);
    dotF(i) = plot3(axF, nan, nan, nan, 'o', 'MarkerSize', 6, 'MarkerFaceColor', c, 'Color', 'k');
    zF(i) = plot3(axF, nan, nan, nan, '-', 'Color', c, 'LineWidth', 2);
end
axis(axF, 'equal'); grid(axF, 'on'); view(axF, 35, 20);
axis(axF, zoomLimits([cell2mat(rF), [pk.viz.collectors, pk.viz.combiner]], 40, 0.2));
set(axF, 'FontSize', 7); xlabel(axF, 'x_f'); ylabel(axF, 'y_f'); zlabel(axF, 'z_f (m)');
title(axF, 'target formation frame O^f (lines: z^i)', 'FontSize', 9, 'FontWeight', 'normal');

% Thruster activity, 5 x 20
axT = axes(fig, 'Position', [0.64 0.37 0.30 0.13]);
img = imagesc(axT, 1:p.nThr, 1:nS, zeros(nS, p.nThr)); hold(axT, 'on');
colormap(axT, flipud(hot)); clim(axT, [0 uMax]);
xline(axT, 5.5:5:15.5, '-', 'Color', [0.5 0.5 0.5]);
set(axT, 'YTick', 1:nS, 'XTick', 3:5:18, 'XTickLabel', compose('pod %d', 1:4), 'FontSize', 7, ...
    'TickLength', [0 0]);
for i = 1:nS, axT.YTickLabel{i} = sprintf('\\color[rgb]{%g,%g,%g}%d', s.sc(i,:), i); end
cb = colorbar(axT); cb.Label.String = 'T/T_{max}'; cb.FontSize = 7;
title(axT, sprintf('thrusters (peak since last frame, scale to fleet peak %.2f)', uMax), ...
      'FontSize', 8, 'FontWeight', 'normal');

% Time histories
axD = axes(fig, 'Position', [0.64 0.235 0.32 0.085]); hold(axD, 'on'); grid(axD, 'on');
axA = axes(fig, 'Position', [0.64 0.14 0.32 0.085]); hold(axA, 'on'); grid(axA, 'on');
axI = axes(fig, 'Position', [0.64 0.045 0.32 0.085]); hold(axI, 'on'); grid(axI, 'on');
patch(axA, [0 tau tau 0], [lo lo hi hi], s.band, 'EdgeColor', 'none');
for i = 1:nS
    r = R{i}; c = s.sc(i,:);
    plot(axD, r.prop.t, max(vecnorm(r.prop.x(8:10,:)), 1e-3), 'Color', c, 'LineWidth', 1.2);
    plot(axA, r.sun.t, r.sun.antiDeg, 'Color', c, 'LineWidth', 1.2);
    plot(axI, r.sol.tNode, r.cumImpulse, 'Color', c, 'LineWidth', 1.2);
end
set(axD, 'YScale', 'log'); ylim(axA, [lo - 5, hi + 5]);
cur = [xline(axD, 0, 'Color', [0.25 0.25 0.25], 'LineWidth', 1.2), ...
       xline(axA, 0, 'Color', [0.25 0.25 0.25], 'LineWidth', 1.2), ...
       xline(axI, 0, 'Color', [0.25 0.25 0.25], 'LineWidth', 1.2)];
xlim([axD axA axI], [0 tau]); set([axD axA], 'XTickLabel', []); set([axD axA axI], 'FontSize', 7);
ylabel(axD, '|r^f_{pi}| (m)'); ylabel(axA, 'anti-Sun (deg)'); ylabel(axI, 'impulse (N s)');
xlabel(axI, 'Time t (s)');

vw = VideoWriter(file, 'MPEG-4'); vw.FrameRate = sched.fps; vw.Quality = sched.quality;
open(vw); cleanup = onCleanup(@() close(vw));
warnState = warning('off', 'MATLAB:hg:DiceyTransformMatrix');
restore = onCleanup(@() warning(warnState));
Lz = 0.08*max(lim(2:2:6) - lim(1:2:5));
for f = 1:numel(sched.t)
    t = sched.t(f);
    C = zeros(nS, p.nThr); Jt = 0;
    for i = 1:nS
        r = R{i};
        ti = min([t, r.prop.t(end), r.sol.tNode(end)]);
        xs = stateAt(r, ti);
        T = thrustInWindow(r, min(sched.t(max(f - 1, 1)), ti), ti);
        C(i,:) = T'/p.Tmax;
        fr = frames(ti, r.p);
        pos = fr.Rbf*(xs(8:10) + r.p.rfp);
        setSpacecraftPose(sc{i}, fr.Rbf*quatToDcm(xs(1:4)), pos, T, plumeLen);
        lbl(i).Position = pos' + [0 0 0.07*(lim(6) - lim(5))];
        upto = r.prop.t <= ti;
        set(trail(i), 'XData', [dB{i}(1,upto) pos(1)], 'YData', [dB{i}(2,upto) pos(2)], ...
            'ZData', [dB{i}(3,upto) pos(3)]);
        tp = fr.Rbf*r.p.rfp;
        set(tgt(i), 'XData', tp(1), 'YData', tp(2), 'ZData', tp(3));
        pf = xs(8:10) + r.p.rfp; z = quatToDcm(xs(1:4))*[0; 0; Lz];
        set(trF(i), 'XData', [rF{i}(1,upto) pf(1)], 'YData', [rF{i}(2,upto) pf(2)], ...
            'ZData', [rF{i}(3,upto) pf(3)]);
        set(dotF(i), 'XData', pf(1), 'YData', pf(2), 'ZData', pf(3));
        set(zF(i), 'XData', pf(1) + [0 z(1)], 'YData', pf(2) + [0 z(2)], 'ZData', pf(3) + [0 z(3)]);
        Jt = Jt + interp1(r.sol.tNode, r.cumImpulse, ti);
    end
    setFormationPose(hf, frames(t, pk).Rbf);
    img.CData = C;
    [cur.Value] = deal(t);
    tr = round(t, 1);
    hTime.String = sprintf('t = %7.1f s (%02d:%02d:%04.1f)   %s   fleet impulse %7.2f N s', t, ...
        floor(tr/3600), floor(mod(tr, 3600)/60), mod(tr, 60), sched.label, Jt);
    writeVideo(vw, grabFrame(fig, sched.height));
end
end
