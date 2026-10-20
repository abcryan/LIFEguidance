function setOrbitLine(h, u, speed, lim)
%SETORBITLINE Draw the orbit line through the origin along unit vector u
%   (scene axes); the label sits inside the axis limits lim (1x6).
u = u/norm(u);
L = 1e4;                                  % far beyond the scene, clipped by the axes
set(h.line, 'XData', [-L L]*u(1), 'YData', [-L L]*u(2), 'ZData', [-L L]*u(3));
s = inf;                                  % largest step along +u inside the box
for i = 1:3
    if u(i) > 1e-9, s = min(s, lim(2*i)/u(i)); elseif u(i) < -1e-9, s = min(s, lim(2*i-1)/u(i)); end
end
q = 0.75*s*u;
set(h.label, 'Position', q', 'String', sprintf('orbital trajectory (%.0f m/s) \\rightarrow', speed));
end
