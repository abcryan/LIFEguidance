function lim = zoomLimits(P, pad, minFrac)
%ZOOMLIMITS Axis limits [xlo xhi ylo yhi zlo zhi] around the points P (3xn)
%   with margin pad; every half-extent is at least minFrac of the largest
%   one, so a nearly planar trajectory still leaves room for the pyramid
%   edges leaving the plane. Objects outside are clipped by the axes.
lo = min(P, [], 2); hi = max(P, [], 2);
c = (lo + hi)/2;
h = (hi - lo)/2 + pad;
h = max(h, minFrac*max(h));
lim = reshape([c - h, c + h]', 1, []);
end
