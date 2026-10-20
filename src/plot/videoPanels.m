function pan = videoPanels(fig, res, p, titleStr)
%VIDEOPANELS Side panels shared by both videos: title, elapsed time and
%   playback rate, twenty-thruster activity bars and a burn timeline.
s = plotStyle();
pan.title = annotation(fig, 'textbox', [0.01 0.935 0.98 0.06], 'String', titleStr, ...
    'EdgeColor', 'none', 'FontSize', 15, 'FontWeight', 'bold', 'HorizontalAlignment', 'left');
pan.time = annotation(fig, 'textbox', [0.66 0.84 0.33 0.09], 'EdgeColor', 'none', ...
    'FontSize', 13, 'FontName', 'FixedWidth', 'FitBoxToText', 'off');
pan.info = annotation(fig, 'textbox', [0.66 0.70 0.33 0.14], 'EdgeColor', 'none', ...
    'FontSize', 10, 'FontName', 'FixedWidth', 'FitBoxToText', 'off', 'Interpreter', 'none');

pan.axT = axes(fig, 'Position', [0.69 0.30 0.29 0.33]);
pan.bar = bar(pan.axT, 1:p.nThr, zeros(1, p.nThr), 'FaceColor', 'flat');
pan.bar.CData = s.pod(p.thr.pod, :);
ylim(pan.axT, [0 1]); xlim(pan.axT, [0.4 p.nThr + 0.6]); grid(pan.axT, 'on');
set(pan.axT, 'XTick', 1:p.nThr, 'XTickLabel', p.thr.label, 'XTickLabelRotation', 90, 'FontSize', 7);
ylabel(pan.axT, 'T_{ik}/T_{max}'); title(pan.axT, 'Thruster activity (colour = pod)', 'FontSize', 10);

pan.axL = axes(fig, 'Position', [0.69 0.10 0.29 0.08]);
hold(pan.axL, 'on');
for b = 1:numel(res.burns)
    B = res.burns(b);
    patch(pan.axL, [B.t0 B.t1 B.t1 B.t0], [0 0 1 1], [0.9 0.2 0.1], 'EdgeColor', [0.9 0.2 0.1]);
end
pan.cursor = xline(pan.axL, 0, 'k', 'LineWidth', 2);
xlim(pan.axL, [0 res.sol.tau]); ylim(pan.axL, [0 1]); set(pan.axL, 'YTick', []);
xlabel(pan.axL, 'Time t (s)'); title(pan.axL, 'Timeline (red: burns)', 'FontSize', 10);
box(pan.axL, 'on');
end
