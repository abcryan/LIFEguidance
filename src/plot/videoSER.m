function videoSER(res, p, sched, file)
%VIDEOSER Animation in the Sun-Earth rotating frame: formation-centred
%   view with SER-parallel axes (rotating virtual formation, its orbital
%   trajectory as a line through the origin, spacecraft path, moving
%   target point, formation axes, plumes). The halo orbit
%   itself is in the static figure halo_orbit.png.
mag = 3; plumeLen = 4;
[xAll, ~] = stateAt(res, res.prop.t);
nP = numel(res.prop.t);
dB = zeros(3, nP); tgtB = zeros(3, nP);
for k = 1:nP
    dB(:,k) = res.Rbf(:,:,k)*(xAll(8:10,k) + p.rfp);
    tgtB(:,k) = res.Rbf(:,:,k)*p.rfp;
end

fig = newFigure('videoSER', [50 50 1280 720]); fig.Resize = 'off';
pan = videoPanels(fig, res, p, ...
    'Sun-Earth rotating frame O^b: formation-centred view, SER-parallel axes');

% Close-up
ax = axes(fig, 'Position', [0.17 0.09 0.49 0.76]); hold(ax, 'on');
hP = plot3(ax, dB(1,:), dB(2,:), dB(3,:), ':', 'Color', [0.4 0.4 0.4]);
hT = plot3(ax, tgtB(1,:), tgtB(2,:), tgtB(3,:), '--', 'Color', [0.2 0.4 0.9]);
trail = plot3(ax, nan, nan, nan, 'k', 'LineWidth', 1.5);
tgtTrail = plot3(ax, nan, nan, nan, 'Color', [0.2 0.4 0.9], 'LineWidth', 2);
tgt = plot3(ax, nan, nan, nan, 'o', 'MarkerSize', 9, 'MarkerFaceColor', [0.2 0.4 0.9], 'Color', 'k');
plot3(ax, 0, 0, 0, 'k+', 'MarkerSize', 10, 'LineWidth', 2);
fAx = gobjects(3, 1); cols = [0.85 0.1 0.1; 0.1 0.6 0.1; 0.1 0.2 0.9];
for i = 1:3
    fAx(i) = plot3(ax, nan, nan, nan, '-', 'Color', cols(i,:), 'LineWidth', 1.2);
end
hf = drawFormation(ax, p);
ho = drawOrbitLine(ax);
h = drawSpacecraft(ax, p, mag);
V = [p.viz.collectors, p.viz.combiner];
VB = zeros(3, 0);
for k = round(linspace(1, nP, 50))
    VB = [VB, res.Rbf(:,:,k)*V]; %#ok<AGROW>
end
lim = sceneLimits([dB, tgtB, VB, zeros(3,1)], 15);
axis(ax, 'equal'); axis(ax, lim); grid(ax, 'on'); view(ax, 35, 20);
camlight(ax, 'headlight'); lighting(ax, 'gouraud');
xlabel(ax, 'x_b (m)'); ylabel(ax, 'y_b (m)'); zlabel(ax, 'z_b (m)');
title(ax, {'origin: formation centre, axes parallel to SER', ...
    sprintf('spacecraft drawn x%d, plume %g m at T_{max}', mag, plumeLen*mag)}, 'FontSize', 9);
legend(ax, [hP, hT, hf.dots, ho.line], {'spacecraft path', 'target path', 'virtual formation', ...
       'orbital trajectory'}, ...
       'Location', 'southoutside', ...
       'Orientation', 'horizontal', 'AutoUpdate', 'off');

vw = VideoWriter(file, 'MPEG-4'); vw.FrameRate = sched.fps; vw.Quality = 75;
open(vw); cleanup = onCleanup(@() close(vw));
warnState = warning('off', 'MATLAB:hg:DiceyTransformMatrix');
restore = onCleanup(@() warning(warnState));
L = 0.25*(lim(2) - lim(1));
for f = 1:numel(sched.t)
    t = sched.t(f);
    x = stateAt(res, t);
    T = thrustInWindow(res, sched.t(max(f - 1, 1)), t);   % peak since last frame
    fr = frames(t, p);
    Rfi = quatToDcm(x(1:4));
    pos = fr.Rbf*(x(8:10) + p.rfp);
    setSpacecraftPose(h, fr.Rbf*Rfi, pos, T, plumeLen);
    upto = res.prop.t <= t;
    set(trail, 'XData', [dB(1,upto) pos(1)], 'YData', [dB(2,upto) pos(2)], 'ZData', [dB(3,upto) pos(3)]);
    tp = fr.Rbf*p.rfp;
    set(tgtTrail, 'XData', [tgtB(1,upto) tp(1)], 'YData', [tgtB(2,upto) tp(2)], 'ZData', [tgtB(3,upto) tp(3)]);
    set(tgt, 'XData', tp(1), 'YData', tp(2), 'ZData', tp(3));
    for i = 1:3
        e = L*fr.Rbf(:,i);
        set(fAx(i), 'XData', [0 e(1)], 'YData', [0 e(2)], 'ZData', [0 e(3)]);
    end
    setFormationPose(hf, fr.Rbf);
    vb = interp1(res.prop.t, res.ref(4:6,:)', t)';          % formation velocity in O^b
    setOrbitLine(ho, vb, norm(vb), lim);
    updatePanels(pan, t, T, sched.label, x, p);
    writeVideo(vw, grabFrame(fig));
end
end

function lim = sceneLimits(P, pad)
lo = min(P, [], 2) - pad; hi = max(P, [], 2) + pad;
lim = reshape([lo hi]', 1, []);
end
