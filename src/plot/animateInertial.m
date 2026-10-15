function animateInertial(log, p, videoFile, opts)
%ANIMATEINERTIAL Video of the maneuver in inertial axes F_e (cf. Fig. 5).
%   The view is centered on the target center of mass with inertial axes.
%   Bodies are drawn to scale; chaser distances beyond opts.compressRadius
%   are compressed logarithmically so the whole approach stays visible.
arguments
    log struct
    p struct
    videoFile (1,:) char
    opts.simStep (1,1) double = 2           % simulated seconds per frame
    opts.frameRate (1,1) double = 30
    opts.compressRadius (1,1) double = 15   % [m], true scale inside
    opts.plumeLength (1,1) double = 3       % [m] at f_max
end

t = log.t;
frames = 1:round(opts.simStep/p.sim.dtLog):numel(t);
rRel = log.s(:,31:33) - log.s(:,13:15);    % r_c/t resolved in F_e
rVis = compressPositions(rRel, opts.compressRadius);

fig = figure('Color', 'w', 'Units', 'pixels', 'Position', [50 50 1280 720], ...
             'Resize', 'off');
ax = axes(fig, 'Position', [0 0.08 0.66 0.84]);
hold(ax, 'on'); axis(ax, 'equal'); axis(ax, 'off');
pad = norm(p.rPt) + 4;
lo = min([min(rVis, [], 1); -pad*[1 1 1]], [], 1) - 3;
hi = max([max(rVis, [], 1);  pad*[1 1 1]], [], 1) + 3;
axis(ax, [lo(1) hi(1) lo(2) hi(2) lo(3) hi(3)]);
view(ax, 35, 18);
axis(ax, 'vis3d'); camzoom(ax, 1.35);
ax.Clipping = 'off';
camlight(ax, 'headlight'); lighting(ax, 'gouraud');
drawInertialTriad(ax, 1.6*norm(p.rPt));

% Target: textured sphere carrying the rock, its LOS cone and body axes
tgt = hgtransform(ax);
drawAsteroid(tgt, norm(p.rPt));
plot3(p.rPt(1), p.rPt(2), p.rPt(3), 'r.', 'MarkerSize', 20, 'Parent', tgt);
drawLosCone(tgt, p, 7);
cols = [0.85 0.1 0.1; 0.1 0.6 0.1; 0.1 0.2 0.9];
for i = 1:3
    e = 1.4*norm(p.rPt)*((1:3)' == i);
    plot3([0 e(1)], [0 e(2)], [0 e(3)], 'Color', cols(i,:), 'LineWidth', 2, ...
          'Parent', tgt);
end

% Chaser: box, end effector, body axes and one plume per thruster
chs = hgtransform(ax);
drawChaserBody(chs, p);
plumes = gobjects(8,1);
[Xc, Yc, Zc] = cylinder([0 0.35], 16);
for k = 1:8
    plumes(k) = hgtransform(chs);
    surf(Xc, Yc, Zc, 'Parent', plumes(k), 'FaceColor', [0.3 0.55 1], ...
         'EdgeColor', 'none', 'FaceAlpha', 0.6);
end
trail = plot3(ax, nan, nan, nan, 'k:', 'LineWidth', 1.3);

hTitle = annotation(fig, 'textbox', [0.02 0.93 0.64 0.06], 'EdgeColor', 'none', ...
                    'FontSize', 14, 'FontWeight', 'bold', ...
                    'HorizontalAlignment', 'center', 'FitBoxToText', 'off');
txt = annotation(fig, 'textbox', [0.02 0.01 0.66 0.06], 'EdgeColor', 'none', ...
                 'FontSize', 11, 'FontName', 'FixedWidth', 'FitBoxToText', 'off');

% Side panels with moving time cursors
axTh = axes(fig, 'Position', [0.70 0.58 0.27 0.32]);
plot(axTh, t, rad2deg(log.theta), 'LineWidth', 1.2); grid(axTh, 'on'); hold(axTh, 'on');
ylabel(axTh, '\theta_{c/t} (deg)'); xlim(axTh, [0 t(end)]);
curTh = plot(axTh, t(1), rad2deg(log.theta(1)), 'ro', 'MarkerFaceColor', 'r');
axRho = axes(fig, 'Position', [0.70 0.13 0.27 0.32]);
plot(axRho, t, log.rho, 'LineWidth', 1.2); grid(axRho, 'on'); hold(axRho, 'on');
ylabel(axRho, '\rho^t_{c/t} (m)'); xlabel(axRho, 'Time (s)'); xlim(axRho, [0 t(end)]);
legend(axRho, {'x', 'y', 'z'}, 'Location', 'northeast', 'AutoUpdate', 'off');
curRho = xline(axRho, 0, 'r', 'LineWidth', 1.2);

vw = VideoWriter(videoFile, 'MPEG-4');
vw.FrameRate = opts.frameRate;
vw.Quality = 95;
open(vw);
cleanup = onCleanup(@() close(vw));
% Thin plume transforms (non-uniform scale) trigger a harmless warning
warnState = warning('off', 'MATLAB:hg:DiceyTransformMatrix');
restoreWarn = onCleanup(@() warning(warnState));

for j = frames
    Rte = projectSO3(reshape(log.s(j,1:9), 3, 3));
    Rce = projectSO3(reshape(log.s(j,19:27), 3, 3));
    tgt.Matrix = blkdiag(Rte, 1);
    chs.Matrix = [Rce rVis(j,:)'; 0 0 0 1];

    F = reshape(log.u(j,:), 3, 8);
    for k = 1:8
        fNorm = norm(F(:,k));
        if fNorm < 1e-2
            plumes(k).Visible = 'off';
            continue
        end
        exhaust = -p.thr.R(:,:,k)*F(:,k)/fNorm;   % opposite to the force
        len = max(opts.plumeLength*fNorm/p.thr.fmax, 0.05);
        plumes(k).Matrix = [alignZ(exhaust) p.thr.pos(:,k); 0 0 0 1] ...
                           * makehgtform('scale', [1 1 len]);
        plumes(k).Visible = 'on';
    end

    set(trail, 'XData', rVis(1:j,1), 'YData', rVis(1:j,2), 'ZData', rVis(1:j,3));
    set(curTh, 'XData', t(j), 'YData', rad2deg(log.theta(j)));
    curRho.Value = t(j);
    hTitle.String = sprintf(['Chaser and target in inertial frame F_e  ' ...
                             '(t = %4.0f s, %gx real time)'], ...
                            t(j), opts.simStep*opts.frameRate);
    txt.String = sprintf(['|rho| = %6.2f m  theta = %6.2f deg  ' ...
                          '|f_c| = %5.1f N  (distances > %g m log-scaled)'], ...
                         norm(log.rho(j,:)), rad2deg(log.theta(j)), ...
                         norm(log.fc(j,:)), opts.compressRadius);
    drawnow;
    writeVideo(vw, getframe(fig));
end
end

function rVis = compressPositions(r, d0)
% Identity inside d0, d0*(1 + log(d/d0)) outside (C1-continuous)
d = vecnorm(r, 2, 2);
scale = ones(size(d));
far = d > d0;
scale(far) = d0*(1 + log(d(far)/d0))./d(far);
rVis = r.*scale;
end

function R = alignZ(v)
% Rotation taking the z-axis to unit vector v
z = [0; 0; 1];
c = z'*v;
if c < -1 + 1e-9
    R = diag([1 -1 -1]);
    return
end
K = skew(cross(z, v));
R = eye(3) + K + K*K/(1 + c);
end

function drawAsteroid(parent, radius)
rng(7);
[X, Y, Z] = sphere(60);
C = zeros(size(X));
for i = 1:25                              % random light/dark surface patches
    c = randn(3,1); c = c/norm(c);
    C = C + sign(randn)*exp(-12*(1 - (c(1)*X + c(2)*Y + c(3)*Z)));
end
surf(radius*X, radius*Y, radius*Z, C, 'Parent', parent, 'EdgeColor', 'none', ...
     'FaceLighting', 'gouraud', 'AmbientStrength', 0.45, ...
     'SpecularStrength', 0.1, 'DiffuseStrength', 0.7);
ax = ancestor(parent, 'axes');
colormap(ax, [linspace(0.15, 0.9, 256)' linspace(0.13, 0.85, 256)' linspace(0.1, 0.75, 256)']);
clim(ax, [-0.5 0.5]*max(abs(C(:))));
% Two body-fixed great circles make the tumbling easy to follow
a = linspace(0, 2*pi, 200);
r = 1.005*radius;
plot3(r*cos(a), r*sin(a), 0*a, 'Color', [0.1 0.1 0.1], 'LineWidth', 1.5, 'Parent', parent);
plot3(r*cos(a), 0*a, r*sin(a), 'Color', [0.1 0.1 0.1], 'LineWidth', 1.5, 'Parent', parent);
end

function drawLosCone(parent, p, height)
tl = tan(p.los.halfAngle);
base = height*[tl tl 1; -tl tl 1; -tl -tl 1; tl -tl 1]' + p.rPt;
for i = 1:4
    j = mod(i, 4) + 1;
    patch([p.rPt(1) base(1,i) base(1,j)], [p.rPt(2) base(2,i) base(2,j)], ...
          [p.rPt(3) base(3,i) base(3,j)], [0.3 0.6 1], 'Parent', parent, ...
          'FaceAlpha', 0.07, 'EdgeColor', [0.3 0.6 1], 'EdgeAlpha', 0.5);
end
end

function drawChaserBody(parent, p)
d = p.dims/2;
V = [-1 -1 -1; 1 -1 -1; 1 1 -1; -1 1 -1; -1 -1 1; 1 -1 1; 1 1 1; -1 1 1].*d';
faces = [1 2 3 4; 5 6 7 8; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8];
patch('Vertices', V, 'Faces', faces, 'FaceColor', [0.2 0.65 0.2], ...
      'EdgeColor', 'k', 'Parent', parent);
plot3([0 p.rPc(1)], [0 p.rPc(2)], [0 p.rPc(3)], 'b-', 'LineWidth', 3, 'Parent', parent);
plot3(p.rPc(1), p.rPc(2), p.rPc(3), 'bo', 'MarkerFaceColor', 'b', ...
      'MarkerSize', 4, 'Parent', parent);
plot3(p.thr.pos(1,:), p.thr.pos(2,:), p.thr.pos(3,:), 'r.', ...
      'MarkerSize', 12, 'Parent', parent);
end

function drawInertialTriad(ax, L)
names = {'x_e', 'y_e', 'z_e'};
for i = 1:3
    e = L*((1:3)' == i);
    plot3(ax, [0 e(1)], [0 e(2)], [0 e(3)], '-', 'Color', [0.6 0.6 0.6]);
    text(ax, 1.05*e(1), 1.05*e(2), 1.05*e(3), names{i}, 'Color', [0.4 0.4 0.4]);
end
end
