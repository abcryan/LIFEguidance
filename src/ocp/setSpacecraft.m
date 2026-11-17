function ocp = setSpacecraft(ocp, pk)
%SETSPACECRAFT The single-spacecraft OCP for the spacecraft of pk (from
%   spacecraftParams): only the NLP parameters, the initial state and the
%   target point, change. The physics the solver was compiled with must
%   match pk (identical spacecraft); otherwise build an OCP from pk.
if isempty(ocp.par)
    error('LIFE:notParametric', 'Build the OCP with buildOcp(p, ref, problem, true).');
end
if ~isequal(physicsKey(pk), ocp.physKey)
    error('LIFE:physicsMismatch', ['Spacecraft %d has other physical parameters than ' ...
          'the compiled OCP; build one with buildOcp(pk, ...).'], pk.sc);
end
ocp.par = [initialRelativeState(pk)./ocp.S; pk.rfp];
end
