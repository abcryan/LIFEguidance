function pan = videoPanels(fig, res, p, titleStr)
%VIDEOPANELS Side panels shared by both videos: title, elapsed time and
%   playback rate, a true-scale body-frame close-up of the spacecraft with
%   its firing nozzles, twenty-thruster activity bars and a burn timeline.
s = plotStyle();
pan.title = annotation(fig, 'textbox', [0.01 0.935 0.98 0.06], 'String', titleStr, ...
    'EdgeColor', 'none', 'FontSize', 15, 'FontWeight', 'bold', 'HorizontalAlignment', 'left');
pan.time = annotation(fig, 'textbox', [0.66 0.86 0.33 0.08], 'EdgeColor', 'none', ...
    'FontSize', 12, 'FontName', 'FixedWidth', 'FitBoxToText', 'off');
pan.info = annotation(fig, 'textbox', [0.66 0.74 0.33 0.12], 'EdgeColor', 'none', ...
    'FontSize', 9, 'FontName', 'FixedWidth', 'FitBoxToText', 'off', 'Interpreter', 'none');

% Body-frame close-up at true scale (fixed camera in O^i)
pan.axS = axes(fig, 'Position', [0.67 0.39 0.32 0.35]); hold(pan.axS, 'on');
pan.sc = drawSpacecraft(pan.axS, p, 1);
axis(pan.axS, 'equal'); axis(pan.axS, [-6.5 6.5 -6.5 6.5 -3 5.5]); view(pan.axS, 35, 25);
cols = [0.85 0.1 0.1; 0.1 0.6 0.1; 0.1 0.2 0.9]; nm = {'x^i', 'y^i', 'z^i'};
for i = 1:3
    e = 6*((1:3)' == i);
    text(pan.axS, e(1), e(2), e(3), nm{i}, 'Color', cols(i,:), 'FontWeight', 'bold', 'FontSize', 9);
end
axis(pan.axS, 'vis3d', 'off'); camzoom(pan.axS, 1.5);
camlight(pan.axS, 'headlight'); lighting(pan.axS, 'gouraud');
annotation(fig, 'textbox', [0.67 0.37 0.32 0.03], 'String', ...
    'spacecraft in body frame O^i, true scale, plume 3 m at T_{max}', 'EdgeColor', 'none', ...
    'FontSize', 8, 'HorizontalAlignment', 'center');

pan.axT = axes(fig, 'Position', [0.70 0.17 0.28 0.17]);
pan.bar = bar(pan.axT, 1:p.nThr, zeros(1, p.nThr), 'FaceColor', 'flat');
pan.bar.CData = s.pod(p.thr.pod, :);
ylim(pan.axT, [0 1]); xlim(pan.axT, [0.4 p.nThr + 0.6]); grid(pan.axT, 'on');
set(pan.axT, 'XTick', 1:p.nThr, 'XTickLabel', p.thr.label, 'XTickLabelRotation', 90, 'FontSize', 7);
ylabel(pan.axT, 'T_{ik}/T_{max}'); title(pan.axT, 'Peak thrust since last frame (colour = pod)', 'FontSize', 9);

pan.axL = axes(fig, 'Position', [0.70 0.03 0.28 0.05]);
hold(pan.axL, 'on');
for b = 1:numel(res.burns)
    B = res.burns(b);
    patch(pan.axL, [B.t0 B.t1 B.t1 B.t0], [0 0 1 1], [0.9 0.2 0.1], 'EdgeColor', [0.9 0.2 0.1]);
end
pan.cursor = xline(pan.axL, 0, 'k', 'LineWidth', 2);
xlim(pan.axL, [0 res.sol.tau]); ylim(pan.axL, [0 1]); set(pan.axL, 'YTick', []);
set(pan.axL, 'FontSize', 7); title(pan.axL, 'Timeline, t (s) (red: burns)', 'FontSize', 9);
box(pan.axL, 'on');
end
