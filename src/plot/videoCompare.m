function videoCompare(R, p, sched, file)
%VIDEOCOMPARE Animation comparing the solutions in the cell array R in the
%   virtual formation frame O^f. The paths nearly coincide (the problems
%   differ in attitude), so each spacecraft gets its own attitude view
%   (top left), seen across the Sun line with anti-Sun up (the camera
%   follows the Sun direction): true-scale body with plumes, the Sun
%   direction, the field of regard of eq. (45) as a green zone of the
%   sphere about the anti-Sun axis, and the body axis z^i, whose tip must
%   stay in it (green inside, red outside). Bottom left: both approach
%   paths. Right: live sky plot of z^i (drawSky), and the anti-Sun angle and
%   cumulative impulse of every problem with a time cursor.
s = plotStyle(); nR = numel(R); plumeLen = 3; rS = 9;   % sphere radius [m]
tMax = max(cellfun(@(r) r.sol.tau, R));
lo = 180 - rad2deg(p.sun.thetaMax); hi = 180 - rad2deg(p.sun.thetaMin);
sunCol = [0.9 0.6 0];

fig = newFigure('videoCompare', [50 50 1280 720]); fig.Resize = 'off';
annotation(fig, 'textbox', [0.01 0.935 0.98 0.06], 'String', ...
    sprintf('Problem %d vs. %d in the virtual formation frame O^f: Sun-angle constraint (45)', ...
    R{1}.problem, R{end}.problem), 'EdgeColor', 'none', 'FontSize', 15, 'FontWeight', 'bold');
hTime = annotation(fig, 'textbox', [0.01 0.885 0.60 0.05], 'EdgeColor', 'none', ...
    'FontSize', 12, 'FontName', 'FixedWidth');

% Attitude views, one per problem
[ph, la] = meshgrid(deg2rad(linspace(lo, hi, 12)), linspace(0, 2*pi, 61));
for i = 1:nR
    c = s.prob(R{i}.problem,:);
    A(i).ax = axes(fig, 'Position', [0.01 + 0.29*(i - 1) 0.42 0.28 0.46]); %#ok<AGROW>
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
    annotation(fig, 'textbox', [0.01 + 0.29*(i - 1) 0.84 0.28 0.04], 'String', ...
        sprintf('Problem %d (true scale)', R{i}.problem), 'Color', c, 'FontSize', 11, ...
        'FontWeight', 'bold', 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
end
annotation(fig, 'textbox', [0.01 0.385 0.57 0.035], 'String', ...
    sprintf(['camera across the Sun line, anti-Sun (dashed) up; green zone: z^i directions ' ...
    '%.0f\\circ-%.0f\\circ from anti-Sun; thick: z^i, green inside, red outside'], lo, hi), 'EdgeColor', 'none', ...
    'FontSize', 8, 'HorizontalAlignment', 'center');

% Approach paths
ax = axes(fig, 'Position', [0.06 0.06 0.50 0.31]); hold(ax, 'on');
rf = cell(1, nR); trail = gobjects(1, nR); dot = gobjects(1, nR);
for i = 1:nR
    c = s.prob(R{i}.problem,:);
    rf{i} = R{i}.prop.x(8:10,:) + p.rfp;
    plot3(ax, rf{i}(1,:), rf{i}(2,:), rf{i}(3,:), ':', 'Color', 0.5 + 0.5*c);
    trail(i) = plot3(ax, nan, nan, nan, '-', 'Color', c, 'LineWidth', 2);
    dot(i) = plot3(ax, nan, nan, nan, 'o', 'MarkerSize', 7, 'MarkerFaceColor', c, 'Color', 'k');
end
plot3(ax, p.rfp(1), p.rfp(2), p.rfp(3), 'o', 'MarkerSize', 9, ...
      'MarkerFaceColor', [0.2 0.4 0.9], 'Color', 'k');
text(ax, p.rfp(1), p.rfp(2), p.rfp(3), '   target point p', 'FontSize', 8, 'Clipping', 'on');
drawFormation(ax, [p.viz.collectors, p.viz.combiner]);
lim = zoomLimits([cell2mat(rf), p.rfp], 15, 0.25);
axis(ax, 'equal'); axis(ax, lim); grid(ax, 'on'); view(ax, 40, 22);
xlabel(ax, 'x_f (m)'); ylabel(ax, 'y_f (m)'); zlabel(ax, 'z_f (m)'); set(ax, 'FontSize', 8);

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
    info(i) = annotation(fig, 'textbox', [0.86 0.76 - 0.17*(i - 1) 0.13 0.13], 'EdgeColor', ...
        s.prob(R{i}.problem,:), 'LineWidth', 1.5, 'FontSize', 9, 'FontName', 'FixedWidth', ...
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
        Rfi = quatToDcm(xs(1:4));
        pos = xs(8:10) + p.rfp;
        % Anti-Sun direction in O^f from the spacecraft (exact, as sunAngle)
        k = min(max(find(r.prop.t <= ti, 1, 'last'), 1), numel(r.prop.t));
        fr = frames(ti, p);
        rbi = r.prop.ref(1:3,k) + fr.Rbf*pos;
        a = fr.Rbf'*(rbi - p.rbs); a = a/norm(a);
        A(i).zone.Matrix = [alignZ(a), zeros(3, 1); 0 0 0 1];
        set(A(i).anti, 'XData', [0 rS*a(1)], 'YData', [0 rS*a(2)], 'ZData', [0 rS*a(3)]);
        set(A(i).sun, 'XData', -[5 11]*a(1), 'YData', -[5 11]*a(2), 'ZData', -[5 11]*a(3));
        A(i).sunTxt.Position = -11.5*a';
        % Sun-pointing camera: anti-Sun up, looking across the Sun line
        v = cross(a, [0; 0; 1]); v = v/norm(v);
        campos(A(i).ax, 40*(cosd(20)*v + sind(20)*a)'); camup(A(i).ax, a');
        camlight(A(i).light, 'headlight');
        setSpacecraftPose(A(i).sc, Rfi, zeros(3, 1), T, plumeLen);
        anti = interp1(r.sun.t, r.sun.antiDeg, ti);
        inBand = anti >= lo - p.gate.sunDeg && anti <= hi + p.gate.sunDeg;
        col = [0.85 0.1 0.1]*~inBand + [0.1 0.7 0.1]*inBand;
        e = rS*Rfi*p.sun.n;
        set(A(i).z, 'XData', [0 e(1)], 'YData', [0 e(2)], 'ZData', [0 e(3)], 'Color', col);
        set(A(i).tip, 'XData', e(1), 'YData', e(2), 'ZData', e(3), 'Color', col);
        upto = r.prop.t <= ti;
        set(trail(i), 'XData', [rf{i}(1,upto) pos(1)], 'YData', [rf{i}(2,upto) pos(2)], ...
            'ZData', [rf{i}(3,upto) pos(3)]);
        set(dot(i), 'XData', pos(1), 'YData', pos(2), 'ZData', pos(3));
        P = interp1(r.sun.t, sky{i}', ti)';
        set(kTrail(i), 'XData', [sky{i}(1,upto) P(1)], 'YData', [sky{i}(2,upto) P(2)]);
        set(kDot(i), 'XData', P(1), 'YData', P(2));
        J = interp1(r.sol.tNode, r.cumImpulse, ti);
        st = 'inside'; if ~inBand, st = 'OUTSIDE'; end
        info(i).String = {sprintf('Problem %d', r.problem), sprintf('anti-Sun %5.1f deg', anti), ...
                          st, sprintf('impulse %6.2f N s', J)};
    end
    curA.Value = t; curI.Value = t;
    tr = round(t, 1);
    hTime.String = sprintf('t = %7.1f s (%02d:%02d:%04.1f)   %s', t, floor(tr/3600), ...
        floor(mod(tr, 3600)/60), mod(tr, 60), sched.label);
    writeVideo(vw, grabFrame(fig, sched.height));
end
end
