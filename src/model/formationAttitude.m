function [R, q] = formationAttitude(spec, t)
%FORMATIONATTITUDE Attitude R^c of a virtual formation frame at time t:
%   skyBasis(spec.sky)*Rz(spec.phase0 + spec.rate*t), and its quaternion
%   (continuous in t). Double or CasADi t (R only for CasADi).
b = spec.phase0 + spec.rate*t;
R = spec.Rstar*rot3(b);
if nargout > 1
    q = quatMul(spec.qstar, [cos(b/2); 0; 0; sin(b/2)]);
end
end
