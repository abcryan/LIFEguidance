function fig = plotHaloOrbit(res, p)
%PLOTHALOORBIT Static view of the reference halo orbit (not animated).
fig = newFigure('halo', [50 50 800 650]);
drawHalo(axes(fig), res, p);
end
