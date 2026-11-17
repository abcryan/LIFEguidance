function [cls, az] = pathClass(varargin)
%PATHCLASS Attitude-path class of a maneuver (Sec. 4.6): the side on which
%   n^i passes the anti-Sun direction. +1 north, -1 south (ecliptic north,
%   as the sky plot; +z_c lies 23 deg from it, on the north side).
%   cls = pathClass(az): from the azimuth of n^i about the Sun line
%   (sunAngle's sun.azimuth, ecliptic north = 0). The class is the sign of
%   cos of the azimuth midway between start and end along the unwrapped
%   path, i.e. where the path crosses the meridian of the anti-Sun axis.
%   [cls, az] = pathClass(t, q, rbf, p): azimuth first, from pointing-frame
%   quaternions q (4xn) at times t and the formation centre rbf (3xn, SER)
%   for the Sun direction (Sec. 5.2), with p in pointing form
%   (pointingParams); used on collocation solutions.
if nargin == 1
    az = varargin{1};
else
    [t, q, rbf, p] = varargin{:};
    n = numel(t); az = zeros(1, n);
    for k = 1:n
        fr = frames(t(k), p);
        a = rbf(:,k) - p.rbs; a = a/norm(a);             % anti-Sun direction
        u1 = [0; 0; 1] - a(3)*a; u1 = u1/norm(u1);       % ecliptic north
        nb = fr.Rbf*quatToDcm(q(:,k))*p.sun.n;
        az(k) = atan2(nb'*cross(a, u1), nb'*u1);
    end
end
azU = unwrap(az);
cls = sign(cos((azU(1) + azU(end))/2));
end
