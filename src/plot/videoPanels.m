function pan = videoPanels(fig, res, p, titleStr)
%VIDEOPANELS Panels shared by both videos: title with elapsed time and
%   playback rate, a true-scale body-frame inset of the spacecraft, and
%   three time histories with a moving cursor (cf. the Malladi animation):
%   eigenaxis attitude error, relative position r^f_pi, and the twenty
%   thrust magnitudes with Tmax.
s = plotStyle();
pan.title = annotation(fig, 'textbox', [0.01 0.935 0.98 0.06], 'String', titleStr, ...
    'EdgeColor', 'none', 'FontSize', 15, 'FontWeight', 'bold', 'HorizontalAlignment', 'left');
pan.time = annotation(fig, 'textbox', [0.01 0.885 0.60 0.05], 'EdgeColor', 'none', ...
    'FontSize', 12, 'FontName', 'FixedWidth', 'FitBoxToText', 'off');

% Body-frame inset at true scale (fixed camera in O^i)
pan.axS = axes(fig, 'Position', [0.00 0.66 0.15 0.23]); hold(pan.axS, 'on');
pan.sc = drawSpacecraft(pan.axS, p, 1);
cols = [0.85 0.1 0.1; 0.1 0.6 0.1; 0.1 0.2 0.9]; nm = {'x^i', 'y^i', 'z^i'};
for i = 1:3
    e = 6*((1:3)' == i);
    text(pan.axS, e(1), e(2), e(3), nm{i}, 'Color', cols(i,:), 'FontWeight', 'bold', 'FontSize', 9);
end
axis(pan.axS, 'equal'); axis(pan.axS, [-6.5 6.5 -6.5 6.5 -3 5.5]); view(pan.axS, 35, 25);
axis(pan.axS, 'vis3d', 'off'); camzoom(pan.axS, 1.4);
camlight(pan.axS, 'headlight'); lighting(pan.axS, 'gouraud');
annotation(fig, 'textbox', [0.00 0.635 0.15 0.03], 'String', ...
    {'body frame O^i, true scale'}, 'EdgeColor', 'none', ...
    'FontSize', 7, 'HorizontalAlignment', 'center');

t = res.prop.t; tau = res.sol.tau;
pan.tHist = t;

% 1) Eigenaxis attitude error with a moving marker
pan.axA = axes(fig, 'Position', [0.70 0.70 0.28 0.19]); hold(pan.axA, 'on'); grid(pan.axA, 'on');
pan.thetaDeg = rad2deg(res.thetaAtan);
plot(pan.axA, t, pan.thetaDeg, 'Color', s.comp(1,:), 'LineWidth', 1.2);
pan.dot = plot(pan.axA, 0, pan.thetaDeg(1), 'o', 'MarkerSize', 7, ...
               'MarkerFaceColor', [0.95 0.2 0.15], 'MarkerEdgeColor', [0.95 0.2 0.15]);
xlim(pan.axA, [0 tau]); ylabel(pan.axA, '\theta_i (deg)'); set(pan.axA, 'FontSize', 8);

% 2) Relative position components with a time cursor
pan.axR = axes(fig, 'Position', [0.70 0.41 0.28 0.19]); hold(pan.axR, 'on'); grid(pan.axR, 'on');
for c = 1:3
    plot(pan.axR, t, res.prop.x(7 + c,:), 'Color', s.comp(c,:), 'LineWidth', 1.2);
end
pan.curR = xline(pan.axR, 0, 'Color', [0.25 0.25 0.25], 'LineWidth', 1.2);
legend(pan.axR, {'x', 'y', 'z'}, 'Location', 'northeast', 'AutoUpdate', 'off', 'FontSize', 7);
xlim(pan.axR, [0 tau]); ylabel(pan.axR, 'r^f_{pi} (m)'); set(pan.axR, 'FontSize', 8);

% 3) All twenty thrust magnitudes (pod colour, five shades per pod) and Tmax
pan.axT = axes(fig, 'Position', [0.70 0.10 0.28 0.21]); hold(pan.axT, 'on'); grid(pan.axT, 'on');
sol = res.sol;
hp = gobjects(1, 4);
for k = 1:p.nThr
    j = p.thr.pod(k); shade = mod(k - 1, 5)/5;
    col = s.pod(j,:)*(1 - 0.6*shade) + 0.6*shade*[1 1 1].*0.8;
    hk = stairs(pan.axT, sol.tNode, [sol.T(k,:), sol.T(k,end)], 'Color', col, 'LineWidth', 1.1);
    if mod(k - 1, 5) == 0, hp(j) = hk; end
end
hmax = yline(pan.axT, p.Tmax, 'k--', 'LineWidth', 1.2);
pan.curT = xline(pan.axT, 0, 'Color', [0.25 0.25 0.25], 'LineWidth', 1.2);
legend(pan.axT, [hp, hmax], {'pod 1', 'pod 2', 'pod 3', 'pod 4', 'T_{max}'}, ...
       'Location', 'east', 'AutoUpdate', 'off', 'FontSize', 7);
xlim(pan.axT, [0 tau]); ylim(pan.axT, [0 1.1*p.Tmax]);
ylabel(pan.axT, 'T_{ik} (N)'); xlabel(pan.axT, 'Time t (s)'); set(pan.axT, 'FontSize', 8);
end
