function videoCompare(R, p, sched, file)
%VIDEOCOMPARE Animation comparing the solutions in the cell array R (each
%   with its own parameters R{i}.p, i.e. its formation phase). Top left:
%   one attitude view per problem in SER axes, seen across the Sun line
%   with anti-Sun up (the same camera for all, following the Sun
%   direction): true-scale body with plumes, the Sun direction, the field
%   of regard of eq. (45) as a green zone of the sphere about the
%   anti-Sun axis, and the body axis z^i, whose tip must stay in it (green
%   inside, red outside). Bottom left: formation-centred close-up with
%   SER-parallel axes: each distinct virtual formation (prescribed phase,
%   Problem 3 phase) with its moving target point, and every path. Right:
%   live sky plot of z^i (drawSky), and the anti-Sun angle and cumulative
%   impulse of every problem with a time cursor.
s = plotStyle(); nR = numel(R); plumeLen = 3; rS = 9;   % sphere radius [m]
tMax = max(cellfun(@(r) r.sol.tau, R));
lo = 180 - rad2deg(p.sun.thetaMax); hi = 180 - rad2deg(p.sun.thetaMin);
sunCol = [0.9 0.6 0];

fig = newFigure('videoCompare', [50 50 1280 720]); fig.Resize = 'off';
annotation(fig, 'textbox', [0.01 0.935 0.98 0.06], 'String', ...
    sprintf('Problems %s: Sun-angle constraint (45) and free arrival phase', ...
    strjoin(cellfun(@(r) sprintf('%d', r.problem), R, 'UniformOutput', false), ', ')), ...
    'EdgeColor', 'none', 'FontSize', 15, 'FontWeight', 'bold');
hTime = annotation(fig, 'textbox', [0.01 0.885 0.60 0.05], 'EdgeColor', 'none', ...
    'FontSize', 12, 'FontName', 'FixedWidth');

% Attitude views, one per problem
[ph, la] = meshgrid(deg2rad(linspace(lo, hi, 12)), linspace(0, 2*pi, 61));
wA = 0.58/nR;
for i = 1:nR
    c = s.prob(R{i}.problem,:);
    A(i).ax = axes(fig, 'Position', [0.01 + wA*(i - 1) 0.42 wA 0.46]); %#ok<AGROW>
    ax = A(i).ax; hold(ax, 'on');
    A(i).zone = hgtransform(ax);
    surf(rS*sin(ph).*cos(la), rS*sin(ph).*sin(la), rS*cos(ph), 'Parent', A(i).zone, ...
         'FaceColor', [0.3 0.8 0.3], 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'FaceLighting', 'none');
    A(i).anti = plot3(ax, nan, nan, nan, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.2);
    A(i).sun = plot3(ax, nan, nan, nan, '-', 'Color', sunCol, 'LineWidth', 3);
    A(i).sunTxt = text(ax, 0, 0, 0, 'Sun', 'Color', sunCol, 'FontWeight', 'bold', 'FontSize', 9);
    A(i).sc = drawSpacecraft(ax, p, 1);
    A(i).z = plot3(ax, nan, nan, nan, '-', 'LineWidth', 4);
    A(i).tip = plot3(ax, nan, nan, nan, 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'w', 'LineWidth', 2);
    axis(ax, 'equal'); axis(ax, 11.5*[-1 1 -1 1 -1 1]); axis(ax, 'off', 'vis3d');
    camva(ax, 36); camtarget(ax, [0 0 0]);
    A(i).light = camlight(ax, 'headlight'); lighting(ax, 'gouraud');
    annotation(fig, 'textbox', [0.01 + wA*(i - 1) 0.84 wA 0.04], 'String', ...
        sprintf('Problem %d (true scale)', R{i}.problem), 'Color', c, 'FontSize', 11, ...
        'FontWeight', 'bold', 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
end
annotation(fig, 'textbox', [0.01 0.385 0.57 0.035], 'String', ...
    sprintf(['SER axes, camera across the Sun line, anti-Sun (dashed) up; green zone: z^i directions ' ...
    '%.0f\\circ-%.0f\\circ from anti-Sun; thick: z^i, green inside, red outside'], lo, hi), 'EdgeColor', 'none', ...
    'FontSize', 8, 'HorizontalAlignment', 'center');

% Formation-centred close-up, SER-parallel axes: one virtual formation per
% distinct phase (Problems 1 and 2 share the prescribed one)
ax = axes(fig, 'Position', [0.06 0.06 0.50 0.30]); hold(ax, 'on');
dB = cell(1, nR); tgtB = cell(1, nR); trail = gobjects(1, nR); dot = gobjects(1, nR);
for i = 1:nR
    c = s.prob(R{i}.problem,:);
    r = R{i}; n = numel(r.prop.t);
    dB{i} = zeros(3, n); tgtB{i} = zeros(3, n);
    for k = 1:n
        dB{i}(:,k) = r.Rbf(:,:,k)*(r.prop.x(8:10,k) + p.rfp);
        tgtB{i}(:,k) = r.Rbf(:,:,k)*p.rfp;
    end
    plot3(ax, dB{i}(1,:), dB{i}(2,:), dB{i}(3,:), ':', 'Color', 0.5 + 0.5*c);
    trail(i) = plot3(ax, nan, nan, nan, '-', 'Color', c, 'LineWidth', 2);
    dot(i) = plot3(ax, nan, nan, nan, 'o', 'MarkerSize', 7, 'MarkerFaceColor', c, 'Color', 'k');
end
phases = cellfun(@(r) r.p.target.phase0, R);
[~, iu] = unique(round(phases, 9), 'stable');
F = struct('h', {}, 'tgt', {}, 'p', {});
hL = gobjects(1, numel(iu)); lblF = cell(1, numel(iu));
for j = 1:numel(iu)
    r = R{iu(j)};
    if abs(r.p.target.phase0 - p.target.phase0) < 1e-9
        col = [0.8 0.1 0.1]; lblF{j} = 'prescribed formation and target point';
    else
        col = s.prob(r.problem,:);
        lblF{j} = sprintf('Problem %d formation (\\phi_\\tau^* = %.1f\\circ)', r.problem, ...
                          rad2deg(angle(exp(1i*r.sol.phiTau))));
    end
    F(j).h = drawFormation(ax, [p.viz.collectors, p.viz.combiner], col);
    F(j).tgt = plot3(ax, nan, nan, nan, 'p', 'MarkerSize', 13, 'MarkerFaceColor', col, 'Color', 'k');
    F(j).p = r.p;
    hL(j) = F(j).tgt;
end
plot3(ax, 0, 0, 0, 'k+', 'MarkerSize', 9, 'LineWidth', 1.5);
lim = zoomLimits([cell2mat(dB), cell2mat(tgtB)], 15, 0.25);
axis(ax, 'equal'); axis(ax, lim); grid(ax, 'on'); view(ax, 35, 20);
xlabel(ax, 'x_b (m)'); ylabel(ax, 'y_b (m)'); zlabel(ax, 'z_b (m)'); set(ax, 'FontSize', 8);
legend(ax, hL, lblF, 'Location', 'northwest', 'FontSize', 7, 'AutoUpdate', 'off');

% Sky plot with full paths (faint) and live traces
axK = axes(fig, 'Position', [0.60 0.42 0.25 0.47]); drawSky(axK, p);
sky = cell(1, nR); kTrail = gobjects(1, nR); kDot = gobjects(1, nR);
for i = 1:nR
    c = s.prob(R{i}.problem,:);
    sky{i} = skyXY(R{i}.sun);
    plot(axK, sky{i}(1,:), sky{i}(2,:), ':', 'Color', 0.5 + 0.5*c);
    kTrail(i) = plot(axK, nan, nan, '-', 'Color', c, 'LineWidth', 1.8);
    kDot(i) = plot(axK, nan, nan, 'o', 'MarkerSize', 8, 'MarkerFaceColor', c, 'Color', 'k');
end
title(axK, 'z^i about the Sun line', 'FontSize', 9, 'FontWeight', 'normal');

% Live numbers per problem
info = gobjects(1, nR);
for i = 1:nR
    info(i) = annotation(fig, 'textbox', [0.855 0.755 - 0.15*(i - 1) 0.14 0.14], 'EdgeColor', ...
        s.prob(R{i}.problem,:), 'LineWidth', 1.5, 'FontSize', 8.5, 'FontName', 'FixedWidth', ...
        'Interpreter', 'none', ...
        'BackgroundColor', 'w', 'Color', s.prob(R{i}.problem,:));
end

% Anti-Sun angle and cumulative impulse with a cursor
axA = axes(fig, 'Position', [0.64 0.22 0.34 0.15]); hold(axA, 'on'); grid(axA, 'on');
sunBands(axA, p, [0 tMax], R{1}.sun.intervals, s);
axI = axes(fig, 'Position', [0.64 0.05 0.34 0.12]); hold(axI, 'on'); grid(axI, 'on');
for i = 1:nR
    c = s.prob(R{i}.problem,:);
    plot(axA, R{i}.sun.t, R{i}.sun.antiDeg, 'Color', c, 'LineWidth', 1.5);
    plot(axI, R{i}.sol.tNode, R{i}.cumImpulse, 'Color', c, 'LineWidth', 1.5);
end
curA = xline(axA, 0, 'Color', [0.25 0.25 0.25], 'LineWidth', 1.2);
curI = xline(axI, 0, 'Color', [0.25 0.25 0.25], 'LineWidth', 1.2);
xlim(axA, [0 tMax]); xlim(axI, [0 tMax]);
ylabel(axA, 'anti-Sun (deg)'); ylabel(axI, 'impulse (N s)'); xlabel(axI, 'Time t (s)');
set([axA axI], 'FontSize', 8); set(axA, 'XTickLabel', []);

vw = VideoWriter(file, 'MPEG-4'); vw.FrameRate = sched.fps; vw.Quality = sched.quality;
open(vw); cleanup = onCleanup(@() close(vw));
warnState = warning('off', 'MATLAB:hg:DiceyTransformMatrix');
restore = onCleanup(@() warning(warnState));
for f = 1:numel(sched.t)
    t = sched.t(f);
    for i = 1:nR
        r = R{i};
        ti = min([t, r.prop.t(end), r.sol.tNode(end)]);   % tau to round-off
        xs = stateAt(r, ti);
        T = thrustInWindow(r, min(sched.t(max(f - 1, 1)), ti), ti);
        if t > r.sol.tau, T = 0*T; end         % arrived: no feedforward thrust
        fr = frames(ti, r.p);
        Rbi = fr.Rbf*quatToDcm(xs(1:4));       % body to SER axes
        pos = fr.Rbf*(xs(8:10) + p.rfp);       % r_bi - r_bf in SER axes
        % Anti-Sun direction from the spacecraft (exact, as sunAngle)
        k = min(max(find(r.prop.t <= ti, 1, 'last'), 1), numel(r.prop.t));
        a = r.prop.ref(1:3,k) + pos - p.rbs; a = a/norm(a);
        A(i).zone.Matrix = [alignZ(a), zeros(3, 1); 0 0 0 1];
        set(A(i).anti, 'XData', [0 rS*a(1)], 'YData', [0 rS*a(2)], 'ZData', [0 rS*a(3)]);
        set(A(i).sun, 'XData', -[5 11]*a(1), 'YData', -[5 11]*a(2), 'ZData', -[5 11]*a(3));
        A(i).sunTxt.Position = -11.5*a';
        % Sun-pointing camera: anti-Sun up, looking across the Sun line
        v = cross(a, [0; 0; 1]); v = v/norm(v);
        campos(A(i).ax, 40*(cosd(20)*v + sind(20)*a)'); camup(A(i).ax, a');
        camlight(A(i).light, 'headlight');
        setSpacecraftPose(A(i).sc, Rbi, zeros(3, 1), T, plumeLen);
        anti = interp1(r.sun.t, r.sun.antiDeg, ti);
        inBand = anti >= lo - p.gate.sunDeg && anti <= hi + p.gate.sunDeg;
        col = [0.85 0.1 0.1]*~inBand + [0.1 0.7 0.1]*inBand;
        e = rS*Rbi*p.sun.n;
        set(A(i).z, 'XData', [0 e(1)], 'YData', [0 e(2)], 'ZData', [0 e(3)], 'Color', col);
        set(A(i).tip, 'XData', e(1), 'YData', e(2), 'ZData', e(3), 'Color', col);
        upto = r.prop.t <= ti;
        set(trail(i), 'XData', [dB{i}(1,upto) pos(1)], 'YData', [dB{i}(2,upto) pos(2)], ...
            'ZData', [dB{i}(3,upto) pos(3)]);
        set(dot(i), 'XData', pos(1), 'YData', pos(2), 'ZData', pos(3));
        P = interp1(r.sun.t, sky{i}', ti)';
        set(kTrail(i), 'XData', [sky{i}(1,upto) P(1)], 'YData', [sky{i}(2,upto) P(2)]);
        set(kDot(i), 'XData', P(1), 'YData', P(2));
        J = interp1(r.sol.tNode, r.cumImpulse, ti);
        st = 'inside'; if ~inBand, st = 'OUTSIDE'; end
        tag = 'prescr.'; if r.problem >= 3, tag = 'free'; end
        info(i).String = {sprintf('Problem %d', r.problem), sprintf('anti-Sun %5.1f deg', anti), ...
                          st, sprintf('impulse %6.2f N s', J), ...
                          sprintf('phase %5.1f (%s)', rad2deg(angle(exp(1i*r.sol.phiTau))), tag)};
    end
    for j = 1:numel(F)                   % formations keep moving after arrival
        frj = frames(t, F(j).p);
        setFormationPose(F(j).h, frj.Rbf);
        tp = frj.Rbf*p.rfp;
        set(F(j).tgt, 'XData', tp(1), 'YData', tp(2), 'ZData', tp(3));
    end
    curA.Value = t; curI.Value = t;
    tr = round(t, 1);
    hTime.String = sprintf('t = %7.1f s (%02d:%02d:%04.1f)   %s', t, floor(tr/3600), ...
        floor(mod(tr, 3600)/60), mod(tr, 60), sched.label);
    writeVideo(vw, grabFrame(fig, sched.height));
end
end
