function P = skyXY(sun, k)
%SKYXY Sky-plot coordinates (2xn) of drawSky for the samples k of sunAngle.
if nargin < 2, k = 1:numel(sun.t); end
P = sun.antiDeg(k).*[sin(sun.azimuth(k)); cos(sun.azimuth(k))];
end
