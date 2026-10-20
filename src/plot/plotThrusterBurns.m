function fig = plotThrusterBurns(res, p)
%PLOTTHRUSTERBURNS Magnified thrust histories around each burn.
sol = res.sol; s = plotStyle();
burns = res.burns; nb = max(numel(burns), 1);
fig = newFigure('burns', [50 50 min(400*nb, 2400) 500]);
tl = tiledlayout(fig, 1, nb, 'TileSpacing', 'compact');
title(tl, 'Thrust magnitudes during each burn (colour = pod)');
for b = 1:numel(burns)
    B = burns(b);
    pad = max(0.25*B.duration, 2);
    ax = nexttile(tl); hold(ax, 'on'); grid(ax, 'on');
    for k = 1:p.nThr
        if max(sol.T(k, B.n0:B.n1)) > p.mesh.burnTol
            stairs(ax, sol.tNode, [sol.T(k,:), sol.T(k,end)], 'Color', s.pod(p.thr.pod(k),:), ...
                   'LineWidth', 1.4, 'DisplayName', char(p.thr.label(k)));
        end
    end
    yline(ax, [0 p.Tmax], 'k--', 'HandleVisibility', 'off');
    xlim(ax, [B.t0 - pad, B.t1 + pad]); ylim(ax, [-0.5 p.Tmax + 0.5]);
    title(ax, sprintf('burn %d: %.1f N s', b, B.impulse), 'FontWeight', 'normal');
    xlabel(ax, 'Time t (s)');
    if b == 1, ylabel(ax, 'T_{ik} (N)'); end
    legend(ax, 'Location', 'best', 'FontSize', 7);
end
end
