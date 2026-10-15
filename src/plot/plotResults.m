function figs = plotResults(log, p, outDir)
%PLOTRESULTS Reproduce Figs. 4-7 of the paper and save them to outDir.
if nargin < 3, outDir = ''; end
t = log.t;
lw = {'LineWidth', 1.2};
xyz = {'$x$', '$y$', '$z$'};

% Fig. 6: attitude/angular-velocity errors, feature point position/velocity
figs.fig6 = figure('Name', 'Fig. 6', 'Color', 'w', 'Position', [100 100 1000 520]);
tl = tiledlayout(2, 2, 'TileSpacing', 'compact');
k = t <= 800;
nexttile; plot(t(k), rad2deg(log.theta(k)), lw{:});
ylabel('$\theta_{c/t}$ (deg)'); grid on;
nexttile; plot(t(k), log.rho(k,:), lw{:});
ylabel('$\rho^t_{c/t}$ (m)'); grid on; legend(xyz, 'Interpreter', 'latex');
nexttile; plot(t(k), log.wct(k,:)*1e3, lw{:});
ylabel('$\omega^t_{c/t}$ ($10^{-3}$ rad/s)'); xlabel('Time (s)'); grid on;
legend(xyz, 'Interpreter', 'latex');
nexttile; plot(t(k), log.rhoDot(k,:), lw{:});
ylabel('$\dot\rho^t_{c/t}$ (m/s)'); xlabel('Time (s)'); grid on;
legend(xyz, 'Interpreter', 'latex');
title(tl, 'Relative attitude and feature-point errors (cf. Fig. 6)');
setLatex(figs.fig6);

% Fig. 7: resultant force/torque and per-thruster magnitude/gimbal angles
figs.fig7 = figure('Name', 'Fig. 7', 'Color', 'w', 'Position', [120 120 1000 620]);
tl = tiledlayout(3, 2, 'TileSpacing', 'compact');
nexttile; stairs(t, log.fc, lw{:}); ylabel('$f_c$ (N)'); grid on;
legend(xyz, 'Interpreter', 'latex');
nexttile; stairs(t, log.fNorm, lw{:}); hold on;
yline(p.thr.fmax, 'k--'); ylabel('$\|f_{t_k}\|$ (N)'); grid on;
nexttile([1 1]); stairs(t, log.tauc, lw{:}); ylabel('$\tau_c$ (N m)'); grid on;
legend(xyz, 'Interpreter', 'latex');
nexttile; stairs(t, log.alpha1, lw{:}); hold on;
yline([-1 1]*p.thr.halfAngle, 'k--'); ylabel('$\alpha_1$ (rad)'); grid on;
k = t >= 800;
nexttile; stairs(t(k), log.fc(k,:), lw{:}); ylabel('$f_c$ (N), zoom');
xlabel('Time (s)'); grid on;
nexttile; stairs(t, log.alpha3, lw{:}); hold on;
yline([-1 1]*p.thr.halfAngle, 'k--'); ylabel('$\alpha_3$ (rad)');
xlabel('Time (s)'); grid on;
title(tl, 'Chaser forces/torques and thruster commands (cf. Fig. 7)');
setLatex(figs.fig7);

% Fig. 4: maneuver in the target frame F_t
figs.fig4 = figure('Name', 'Fig. 4', 'Color', 'w', 'Position', [140 140 800 650]);
ax = axes; hold(ax, 'on');
rct = zeros(numel(t), 3);
for j = 1:numel(t)
    Rct = reshape(log.z(j,1:9), 3, 3);
    rct(j,:) = (log.rho(j,:)' - Rct*p.rPc + p.rPt)';
end
drawAsteroid(ax, eye(3), zeros(3,1), norm(p.rPt));
plot3(ax, p.rPt(1), p.rPt(2), p.rPt(3), 'r.', 'MarkerSize', 18);
drawLosCone(ax, p, 60);
plot3(ax, rct(:,1), rct(:,2), rct(:,3), 'k:', 'LineWidth', 1.2);
for j = 1:120:numel(t)
    drawChaser(ax, p, reshape(log.z(j,1:9), 3, 3), rct(j,:)', 1);
end
drawChaser(ax, p, reshape(log.z(end,1:9), 3, 3), rct(end,:)', 1);
axis(ax, 'equal'); grid(ax, 'on'); view(ax, 120, 20);
xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
title('Chaser and target in target frame F_t (cf. Fig. 4)');

% Fig. 5: same maneuver seen in inertial axes (centered on the target)
figs.fig5 = figure('Name', 'Fig. 5', 'Color', 'w', 'Position', [160 160 800 650]);
ax = axes; hold(ax, 'on');
scale = 0.1;                            % "scaled for visibility"
rRel = log.s(:,31:33) - log.s(:,13:15);
plot3(ax, scale*rRel(:,1), scale*rRel(:,2), scale*rRel(:,3), 'k:', 'LineWidth', 1.2);
drawAsteroid(ax, eye(3), zeros(3,1), 0.5*norm(p.rPt)*scale*4);
for j = [find(t <= 600 & mod(t, 120) == 0)' numel(t)]
    Rte = reshape(log.s(j,1:9), 3, 3);
    Rce = reshape(log.s(j,19:27), 3, 3);
    pos = scale*rRel(j,:)';
    drawChaser(ax, p, Rce, pos, 0.25);
    drawTargetTriad(ax, Rte, pos);
end
axis(ax, 'equal'); grid(ax, 'on'); view(ax, 35, 25);
xlabel('x'); ylabel('y'); zlabel('z');
title({'Chaser in inertial axes F_e, relative to target (cf. Fig. 5)', ...
       'positions scaled by 0.1, chaser by 0.25; triads show target attitude'});

if ~isempty(outDir)
    if ~isfolder(outDir), mkdir(outDir); end
    names = fieldnames(figs);
    for i = 1:numel(names)
        exportgraphics(figs.(names{i}), fullfile(outDir, [names{i} '.png']), ...
                       'Resolution', 150);
    end
end
end

function setLatex(fig)
set(findall(fig, '-property', 'Interpreter'), 'Interpreter', 'latex');
set(findall(fig, '-property', 'TickLabelInterpreter'), 'TickLabelInterpreter', 'latex');
end

function drawAsteroid(ax, R, c, radius)
[X, Y, Z] = sphere(30);
P = R*[X(:) Y(:) Z(:)]'*radius + c;
surf(ax, reshape(P(1,:), size(X)), reshape(P(2,:), size(X)), reshape(P(3,:), size(X)), ...
     'FaceColor', [0.6 0.6 0.6], 'EdgeColor', 'none', 'FaceAlpha', 0.8);
light(ax, 'Position', [1 -1 2]); lighting(ax, 'gouraud');
end

function drawLosCone(ax, p, height)
tl = tan(p.los.halfAngle);
base = height*[tl tl 1; -tl tl 1; -tl -tl 1; tl -tl 1]' + p.rPt;
for i = 1:4
    j = mod(i, 4) + 1;
    fill3(ax, [p.rPt(1) base(1,i) base(1,j)], [p.rPt(2) base(2,i) base(2,j)], ...
          [p.rPt(3) base(3,i) base(3,j)], [0.3 0.6 1], ...
          'FaceAlpha', 0.08, 'EdgeColor', [0.3 0.6 1]);
end
end

function drawChaser(ax, p, R, c, scale)
d = scale*p.dims/2;
V = ([-1 -1 -1; 1 -1 -1; 1 1 -1; -1 1 -1; -1 -1 1; 1 -1 1; 1 1 1; -1 1 1].*d')';
V = (R*V + c)';
faces = [1 2 3 4; 5 6 7 8; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8];
patch(ax, 'Vertices', V, 'Faces', faces, 'FaceColor', [0.2 0.7 0.2], ...
      'EdgeColor', 'k', 'FaceAlpha', 0.9);
ee = c + R*p.rPc*scale;
plot3(ax, [c(1) ee(1)], [c(2) ee(2)], [c(3) ee(3)], 'b-', 'LineWidth', 2);
end

function drawTargetTriad(ax, R, c)
L = 2;
cols = 'rgb';
for i = 1:3
    e = c + L*R(:,i);
    plot3(ax, [c(1) e(1)], [c(2) e(2)], [c(3) e(3)], cols(i), 'LineWidth', 1);
end
end
