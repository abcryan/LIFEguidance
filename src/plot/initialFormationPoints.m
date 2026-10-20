function V = initialFormationPoints(p)
%INITIALFORMATIONPOINTS The five points of the initial virtual formation
%   at t = 0, resolved in the target formation frame O^f (3x5) [m]; the
%   two frames share their origin at t = 0.
Rf_f0 = formationAttitude(p.target, 0)'*formationAttitude(p.init, 0);
V = Rf_f0*pyramidPoints(p.init);
end
