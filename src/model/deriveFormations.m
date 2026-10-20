function p = deriveFormations(p)
%DERIVEFORMATIONS Quantities derived from the formation specs p.target and
%   p.init: normalized sky directions, skyBasis rotations, the target spin
%   p.ws/p.phi0, the target point p.rfp and the drawn pyramid p.viz.
%   Call again after editing p.target or p.init.
for name = {'target', 'init'}
    spec = p.(name{1});
    spec.sky = spec.sky(:)/norm(spec.sky);
    spec.Rstar = skyBasis(spec.sky);
    spec.qstar = dcmToQuat(spec.Rstar);
    p.(name{1}) = spec;
end
p.ws   = p.target.rate;                 % target formation spin [rad/s]
p.phi0 = p.target.phase0;
P = pyramidPoints(p.target);
p.rfp = P(:, p.target.point);           % target point in O^f [m]
p.viz.collectors = P(:, 1:4);           % drawn in the videos
p.viz.combiner   = P(:, 5);
end
