function pan = videoPanels(fig, res, p, titleStr)
%VIDEOPANELS Panels shared by both videos: title with elapsed time and
%   playback rate, a true-scale body-frame inset of the spacecraft, and
%   three time histories with a moving cursor (cf. the Malladi animation):
%   eigenaxis attitude error, relative position r^f_pi, and the thrust
%   magnitude of each of the twenty thrusters with Tmax.
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
pan.axA = axes(fig, 'Position', [0.70 0.73 0.28 0.17]); hold(pan.axA, 'on'); grid(pan.axA, 'on');
pan.thetaDeg = rad2deg(res.thetaAtan);
plot(pan.axA, t, pan.thetaDeg, 'Color', s.comp(1,:), 'LineWidth', 1.2);
pan.dot = plot(pan.axA, 0, pan.thetaDeg(1), 'o', 'MarkerSize', 7, ...
               'MarkerFaceColor', [0.95 0.2 0.15], 'MarkerEdgeColor', [0.95 0.2 0.15]);
xlim(pan.axA, [0 tau]); ylabel(pan.axA, '\theta_i (deg)'); set(pan.axA, 'FontSize', 8);

% 2) Relative position components with a time cursor
pan.axR = axes(fig, 'Position', [0.70 0.48 0.28 0.17]); hold(pan.axR, 'on'); grid(pan.axR, 'on');
for c = 1:3
    plot(pan.axR, t, res.prop.x(7 + c,:), 'Color', s.comp(c,:), 'LineWidth', 1.2);
end
pan.curR = xline(pan.axR, 0, 'Color', [0.25 0.25 0.25], 'LineWidth', 1.2);
legend(pan.axR, {'x', 'y', 'z'}, 'Location', 'northeast', 'AutoUpdate', 'off', 'FontSize', 7);
xlim(pan.axR, [0 tau]); ylabel(pan.axR, 'r^f_{pi} (m)'); set(pan.axR, 'FontSize', 8);

% 3) Thrust magnitude of every thruster as one vertical line per control
%    interval (at its midpoint), overlaid; the thrusters with the largest
%    total impulse are drawn first so weaker ones stay visible on top
pan.axT = axes(fig, 'Position', [0.70 0.08 0.28 0.20]); hold(pan.axT, 'on'); grid(pan.axT, 'on');
sol = res.sol;
tMid = (sol.tNode(1:end-1) + sol.tNode(2:end))/2;
[~, order] = sort(res.thrusterImpulse, 'descend');
hk = gobjects(1, p.nThr);
for k = order(:)'
    on = sol.T(k,:) > p.burnTol;
    X = [tMid(on); tMid(on); nan(1, nnz(on))];
    Y = [zeros(1, nnz(on)); sol.T(k, on); nan(1, nnz(on))];
    if ~any(on), X = nan; Y = nan; end   % never fires: legend entry only
    hk(k) = plot(pan.axT, X(:), Y(:), '-', 'Color', s.thr(k,:), 'LineWidth', 1.5);
end
hmax = yline(pan.axT, p.Tmax, 'k--', 'LineWidth', 1.2);
pan.curT = xline(pan.axT, 0, 'Color', [0.25 0.25 0.25], 'LineWidth', 1.2);
lg = legend(pan.axT, [hk, hmax], [cellstr(p.thr.label), {'T_{max}'}], 'NumColumns', 5, ...
       'AutoUpdate', 'off', 'FontSize', 6, 'Box', 'off');
lg.ItemTokenSize = [10 6];
lg.Position = [0.69 0.30 0.30 0.11];
xlim(pan.axT, [0 tau]); ylim(pan.axT, [0 1.1*p.Tmax]);
ylabel(pan.axT, 'T_{ik} (N)'); xlabel(pan.axT, 'Time t (s)'); set(pan.axT, 'FontSize', 8);
end
