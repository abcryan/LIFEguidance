function videoCollision(F, P4, p, sched, file)
%VIDEOCOLLISION Problem 5 animation of the collision avoidance (Sec. 5.3),
%   next to Problem 4, same camera, scale and time:
%   top     the four collectors in the pointing frame O* (the combiner stays
%           more than 100 m away and is left out), Problem 4 left, Problem
%           5 right. Every collector carries a translucent sphere of radius
%           dmin/2 at the scale of the trajectories (two spheres touch
%           exactly at d_ij = dmin) and its z^i axis. Problem 4: spheres of
%           a pair closer than dmin turn red, joined by a "VIOLATION" line.
%           Problem 5: a pair at the imposed bound (d_ij < dmin + margin +
%           0.1 m) is joined by a red line, pairs within dmin + nearLabel by
%           an orange one.
%   middle  the anti-Sun angle of z^i of all five spacecraft per problem in
%           the field of regard (Sun-angle constraint (48)).
%   bottom  min over all pairs of d_ij(t), both problems, dmin, cursor.
%   The two problems arrive at slightly different points, since their
%   optimal arrival phases differ.
s = plotStyle(); cd = p.collision; nC = 4; nS = numel(F.res);
R = {P4.fleet.res, F.res}; tau = F.res{1}.sol.tau;
J = [P4.fleet.info.J, F.info.J]; phi = rad2deg([P4.fleet.info.phiTau, F.info.phiTau]);
pos = cell(2, nC); zax = cell(2, nC); tt = cell(1, 2);
for m = 1:2
    tt{m} = R{m}{1}.prop.t;
    for i = 1:nC
        xP = R{m}{i}.prop.xP;
        pos{m, i} = xP(8:10, :);
        z = zeros(3, size(xP, 2));
        for k = 1:size(xP, 2), Rq = quatToDcm(xP(1:4, k)); z(:, k) = Rq(:, 3); end
        zax{m, i} = z;
    end
end
at = @(P, tq, t) interp1(tq, P', min(t, tq(end)))';
pairs = nchoosek(1:nC, 2); nP = size(pairs, 1);
Dmin = cell(1, 2);
for m = 1:2
    D = zeros(nP, numel(tt{m}));
    for k = 1:nP, D(k, :) = vecnorm(pos{m, pairs(k,1)} - pos{m, pairs(k,2)}); end
    Dmin{m} = min(D, [], 1);
end
lo = 180 - rad2deg(p.sun.thetaMax); hi = 180 - rad2deg(p.sun.thetaMin);

fig = newFigure('videoCollision', [50 50 1280 720]); fig.Resize = 'off';
annotation(fig, 'textbox', [0.01 0.94 0.98 0.055], 'String', sprintf(['Collision avoidance: Problem 4 ' ...
    '(no constraint) vs Problem 5 (d_{min} = %.0f m between all spacecraft), \\Pi_5 = %+.2f N s'], ...
    cd.dmin, J(2) - J(1)), 'EdgeColor', 'none', 'FontSize', 15, 'FontWeight', 'bold');
hTime = annotation(fig, 'textbox', [0.01 0.895 0.38 0.045], 'EdgeColor', 'none', 'FontSize', 12, ...
                   'FontName', 'FixedWidth');

allP = cell2mat(pos(:)');
lim = [min(allP, [], 2) - 15, max(allP, [], 2) + 15]';
[xs, ys, zs] = sphere(24); rs = cd.dmin/2; Lz = 2.2*rs;
near = cd.dmin + cd.nearLabel; activeD = cd.dmin + cd.margin + 0.1;
ttl = {sprintf('Problem 4: J_4^* = %.2f N s, \\phi_\\tau^* = %.2f\\circ, min. separation %.2f m', ...
               J(1), phi(1), min(Dmin{1})), ...
       sprintf('Problem 5: J_5^* = %.2f N s, \\phi_\\tau^* = %.2f\\circ, min. separation %.2f m', ...
               J(2), phi(2), min(Dmin{2}))};
ax = gobjects(1, 2); sph = gobjects(2, nC); trail = gobjects(2, nC); lbl = gobjects(2, nC);
zl = gobjects(2, nC); link = gobjects(2, nP); linkTxt = gobjects(2, nP); curA = gobjects(1, 2);
for m = 1:2
    ax(m) = axes(fig, 'Position', [0.025 + 0.49*(m - 1), 0.44, 0.46, 0.39]); hold(ax(m), 'on');
    for i = 1:nC
        c = s.sc(i,:);
        plot3(ax(m), pos{m,i}(1,:), pos{m,i}(2,:), pos{m,i}(3,:), ':', 'Color', 0.4 + 0.6*c);
        trail(m,i) = plot3(ax(m), nan, nan, nan, '-', 'Color', c, 'LineWidth', 2);
        sph(m,i) = surf(ax(m), xs, ys, zs, 'FaceColor', c, 'FaceAlpha', 0.35, 'EdgeColor', 'none');
        zl(m,i) = plot3(ax(m), nan, nan, nan, '-', 'Color', 0.6*c, 'LineWidth', 2);
        lbl(m,i) = text(ax(m), 0, 0, 0, sprintf('%d', i), 'Color', c, 'FontWeight', 'bold', ...
                        'FontSize', 11, 'HorizontalAlignment', 'center');
    end
    for k = 1:nP
        link(m,k) = plot3(ax(m), nan, nan, nan, '-', 'LineWidth', 2.5);
        linkTxt(m,k) = text(ax(m), 0, 0, 0, '', 'FontWeight', 'bold', 'FontSize', 10, ...
                            'BackgroundColor', 'w', 'Margin', 1, 'HorizontalAlignment', 'center');
    end
    axis(ax(m), 'equal'); axis(ax(m), lim(:)'); grid(ax(m), 'on'); view(ax(m), 30, 25);
    axis(ax(m), 'vis3d'); camzoom(ax(m), 1.25);
    camlight(ax(m), 'headlight'); lighting(ax(m), 'gouraud');
    xlabel(ax(m), 'x_* (m)'); ylabel(ax(m), 'y_* (m)'); zlabel(ax(m), 'z_* (m)');
    set(ax(m), 'FontSize', 8);
    annotation(fig, 'textbox', [0.025 + 0.49*(m - 1), 0.85, 0.46, 0.04], 'String', ttl{m}, ...
        'EdgeColor', 'none', 'FontSize', 10.5, 'HorizontalAlignment', 'center');

    % Sun-angle constraint of all five spacecraft
    axA = axes(fig, 'Position', [0.06 + 0.49*(m - 1), 0.255, 0.41, 0.085]); hold(axA, 'on'); grid(axA, 'on');
    patch(axA, [0 tau tau 0], [lo lo hi hi], s.band, 'EdgeColor', 'none');
    for i = 1:nS
        sn = R{m}{i}.sun;
        plot(axA, sn.t, sn.antiDeg, '-', 'Color', s.sc(i,:), 'LineWidth', 1.2);
    end
    curA(m) = xline(axA, 0, 'Color', [0.2 0.2 0.8], 'LineWidth', 1.2);
    xlim(axA, [0 tau]); ylim(axA, [lo - 8, hi + 8]); set(axA, 'FontSize', 7, 'XTickLabel', []);
    ylabel(axA, 'anti-Sun (deg)');
    maxV = max(cellfun(@(r) r.sun.maxViolDeg, R{m}));
    title(axA, sprintf('Sun angle of z^i, all 5 spacecraft (green: field of regard %g-%g\\circ, max. violation %.4f\\circ)', ...
          lo, hi, maxV), 'FontSize', 8, 'FontWeight', 'normal');
end
annotation(fig, 'textbox', [0.40 0.893 0.59 0.045], 'EdgeColor', 'none', 'FontSize', 9, ...
    'HorizontalAlignment', 'right', 'String', sprintf(['collectors in O^*, same scale and camera; spheres ' ...
    'of radius d_{min}/2 = %.0f m to scale (touch at d_{min}); short lines: z^i; combiner not shown'], rs));

% Minimum separation over time
axM = axes(fig, 'Position', [0.06 0.065 0.9 0.14]); hold(axM, 'on'); grid(axM, 'on');
plot(axM, tt{1}, Dmin{1}, '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 1.6);
plot(axM, tt{2}, Dmin{2}, '-', 'Color', 'k', 'LineWidth', 1.8);
yline(axM, cd.dmin, 'r-', 'LineWidth', 1.4);
cur = xline(axM, 0, 'Color', [0.2 0.2 0.8], 'LineWidth', 1.4);
xlim(axM, [0 tau]); ylim(axM, [0 1.05*max([Dmin{:}])]); set(axM, 'FontSize', 8);
xlabel(axM, 'Time t (s)'); ylabel(axM, 'min_{ij} d_{ij} (m)');
legend(axM, {'Problem 4', 'Problem 5', sprintf('d_{min} = %.0f m', cd.dmin)}, 'Location', 'east', ...
       'FontSize', 8, 'AutoUpdate', 'off');

vw = VideoWriter(file, 'MPEG-4'); vw.FrameRate = sched.fps; vw.Quality = sched.quality;
open(vw); cleanup = onCleanup(@() close(vw));
red = [0.85 0.05 0.05]; orange = [0.95 0.55 0.1];
for f = 1:numel(sched.t)
    t = sched.t(f);
    for m = 1:2
        q = cell(1, nC); hot = false(1, nC);
        upto = tt{m} <= t;
        for i = 1:nC
            q{i} = at(pos{m,i}, tt{m}, t); z = at(zax{m,i}, tt{m}, t); z = z/norm(z);
            set(sph(m,i), 'XData', q{i}(1) + rs*xs, 'YData', q{i}(2) + rs*ys, 'ZData', q{i}(3) + rs*zs);
            set(trail(m,i), 'XData', [pos{m,i}(1,upto) q{i}(1)], 'YData', [pos{m,i}(2,upto) q{i}(2)], ...
                'ZData', [pos{m,i}(3,upto) q{i}(3)]);
            set(zl(m,i), 'XData', q{i}(1) + [0 Lz*z(1)], 'YData', q{i}(2) + [0 Lz*z(2)], ...
                'ZData', q{i}(3) + [0 Lz*z(3)]);
            lbl(m,i).Position = q{i}' + [0 0 -2.2*rs];
        end
        for k = 1:nP
            a = q{pairs(k,1)}; b = q{pairs(k,2)}; dk = norm(a - b);
            if dk < near
                if m == 1 && dk < cd.dmin
                    col = red; str = sprintf('VIOLATION %.1f m', dk); hot(pairs(k,:)) = true;
                elseif m == 2 && dk < activeD
                    col = red; str = sprintf('%.1f m (at bound)', dk);
                else
                    col = orange; str = sprintf('%.1f m', dk);
                end
                set(link(m,k), 'XData', [a(1) b(1)], 'YData', [a(2) b(2)], 'ZData', [a(3) b(3)], 'Color', col);
                set(linkTxt(m,k), 'Position', (a' + b')/2 + [0 0 3*rs], 'String', str, ...
                    'Color', col, 'EdgeColor', col);
            else
                set(link(m,k), 'XData', nan, 'YData', nan, 'ZData', nan);
                set(linkTxt(m,k), 'String', '');
            end
        end
        for i = 1:nC                    % overlapping spheres (Problem 4) in red
            if hot(i), set(sph(m,i), 'FaceColor', red, 'FaceAlpha', 0.55);
            else, set(sph(m,i), 'FaceColor', s.sc(i,:), 'FaceAlpha', 0.35); end
        end
        curA(m).Value = t;
    end
    cur.Value = t;
    tr = round(t, 1);
    hTime.String = sprintf('t = %7.1f s (%02d:%02d:%04.1f)   %s', t, floor(tr/3600), ...
        floor(mod(tr, 3600)/60), mod(tr, 60), sched.label);
    writeVideo(vw, grabFrame(fig, sched.height));
end
end
