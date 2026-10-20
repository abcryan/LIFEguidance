function videoSER(res, p, sched, file)
%VIDEOSER Animation in the Sun-Earth rotating frame: halo orbit with the
%   arc flown, and a formation-centred close-up with SER-parallel axes
%   (spacecraft path, moving target point, formation axes, plumes).
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
    'Sun-Earth rotating frame O^b: halo orbit (left) and formation-centred close-up, SER-parallel axes');

% Halo overview, relative to Earth [10^3 km]
axH = axes(fig, 'Position', [0.03 0.52 0.20 0.33]); hold(axH, 'on');
rel = @(r) (r - p.rbe/p.dse)*p.dse/1e6;
H = rel(res.haloFull(1:3,:));
A = rel(res.ref(1:3,:)/p.dse);
plot3(axH, H(1,:), H(2,:), H(3,:), 'Color', [0.6 0.6 0.6], 'LineWidth', 0.5);
plot3(axH, A(1,:), A(2,:), A(3,:), 'r', 'LineWidth', 3);
mk = plot3(axH, A(1,1), A(2,1), A(3,1), 'ko', 'MarkerFaceColor', 'y', 'MarkerSize', 7);
plot3(axH, 0, 0, 0, 'bo', 'MarkerFaceColor', 'b');
text(axH, 0, 0, 0, '  Earth', 'FontSize', 8);
axis(axH, 'equal'); grid(axH, 'on'); view(axH, 35, 25);
xlabel(axH, 'x_b (10^3 km)'); ylabel(axH, 'y_b (10^3 km)'); zlabel(axH, 'z_b (10^3 km)');
title(axH, 'L_2 halo, red: 4 h arc', 'FontSize', 9);
set(axH, 'FontSize', 7);

% Close-up
ax = axes(fig, 'Position', [0.30 0.08 0.36 0.78]); hold(ax, 'on');
plot3(ax, dB(1,:), dB(2,:), dB(3,:), ':', 'Color', [0.4 0.4 0.4]);
plot3(ax, tgtB(1,:), tgtB(2,:), tgtB(3,:), '--', 'Color', [0.2 0.4 0.9]);
trail = plot3(ax, nan, nan, nan, 'k', 'LineWidth', 1.5);
tgtTrail = plot3(ax, nan, nan, nan, 'Color', [0.2 0.4 0.9], 'LineWidth', 2);
tgt = plot3(ax, nan, nan, nan, 'o', 'MarkerSize', 9, 'MarkerFaceColor', [0.2 0.4 0.9], 'Color', 'k');
plot3(ax, 0, 0, 0, 'k+', 'MarkerSize', 10, 'LineWidth', 2);
fAx = gobjects(3, 1); cols = [0.85 0.1 0.1; 0.1 0.6 0.1; 0.1 0.2 0.9];
for i = 1:3
    fAx(i) = plot3(ax, nan, nan, nan, '-', 'Color', cols(i,:), 'LineWidth', 1.2);
end
h = drawSpacecraft(ax, p, mag);
lim = sceneLimits([dB, tgtB, zeros(3,1)], 15);
axis(ax, 'equal'); axis(ax, lim); grid(ax, 'on'); view(ax, 35, 20);
camlight(ax, 'headlight'); lighting(ax, 'gouraud');
xlabel(ax, 'x_b (m)'); ylabel(ax, 'y_b (m)'); zlabel(ax, 'z_b (m)');
title(ax, {'origin: formation centre, axes parallel to SER', ...
    sprintf('spacecraft drawn x%d, plume %g m at T_{max}', mag, plumeLen*mag)}, 'FontSize', 9);
legend(ax, {'spacecraft path', 'target path'}, 'Location', 'southoutside', ...
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
    a = rel(interp1(res.prop.t, res.ref(1:3,:)', t)'/p.dse);
    set(mk, 'XData', a(1), 'YData', a(2), 'ZData', a(3));
    updatePanels(pan, t, T, sched.label, x, p);
    writeVideo(vw, grabFrame(fig));
end
end

function lim = sceneLimits(P, pad)
lo = min(P, [], 2) - pad; hi = max(P, [], 2) + pad;
lim = reshape([lo hi]', 1, []);
end
