function drawSky(ax, p)
%DRAWSKY Sky plot of the direction of n^i about the Sun line: radius =
%   anti-Sun angle (deg, anti-Sun at the centre, Sun on the outer circle),
%   azimuth from ecliptic north (up). The field of regard of eq. (48) is
%   the shaded annulus. Map directions with skyXY.
s = plotStyle();
hold(ax, 'on');
lo = 180 - rad2deg(p.sun.thetaMax); hi = 180 - rad2deg(p.sun.thetaMin);
th = linspace(0, 2*pi, 181);
patch(ax, [hi*cos(th), lo*cos(fliplr(th))], [hi*sin(th), lo*sin(fliplr(th))], s.band, ...
      'EdgeColor', 'none', 'HandleVisibility', 'off');
for r = 30:30:180
    plot(ax, r*cos(th), r*sin(th), '-', 'Color', [0.8 0.8 0.8], 'HandleVisibility', 'off');
    text(ax, r*cos(pi/3), r*sin(pi/3), sprintf('%d\\circ', r), 'FontSize', 7, 'Color', [0.5 0.5 0.5]);
end
plot(ax, [lo*cos(th); hi*cos(th)]', [lo*sin(th); hi*sin(th)]', 'g-', 'HandleVisibility', 'off');
plot(ax, [-180 180], [0 0], ':', [0 0], [-180 180], ':', 'Color', [0.7 0.7 0.7], 'HandleVisibility', 'off');
plot(ax, 0, 0, 'p', 'MarkerSize', 10, 'MarkerFaceColor', 'k', 'Color', 'k', 'HandleVisibility', 'off');
text(ax, 4, -8, 'anti-Sun', 'FontSize', 7);
text(ax, 0, 186, 'ecliptic north', 'FontSize', 7, 'HorizontalAlignment', 'center');
text(ax, 128, -128, 'Sun (180\circ)', 'FontSize', 7, 'Color', [0.8 0.5 0]);
axis(ax, 'equal'); axis(ax, [-195 195 -195 195]); axis(ax, 'off');
end
