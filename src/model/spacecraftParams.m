function pk = spacecraftParams(p, i)
%SPACECRAFTPARAMS Parameters of spacecraft i of the fleet (Problem 4):
%   p with the initial and target point set to p.fleet.point(i) (fixed
%   assignment, the same point of both pyramids), so initialRelativeState,
%   targetState, p.rfp and all post-processing refer to spacecraft i.
%   Spacecraft-specific values in p.fleet.override{i} replace those of p
%   (mass, inertia and SRP area enter Jinv and beta, which are recomputed).
pk = p;
pk.init.point = p.fleet.point(i);
pk.target.point = p.fleet.point(i);
ov = p.fleet.override{i};
if ~isempty(ov)
    for f = fieldnames(ov)'
        pk.(f{1}) = ov.(f{1});
    end
    pk.Jinv = inv(pk.J);
    pk.beta = pk.Ssun*pk.Cr*pk.Asrp*pk.dse^2/(pk.clight*pk.m);
end
pk = deriveFormations(pk);
pk.sc = i;
pk.scName = p.fleet.name(i);
end
