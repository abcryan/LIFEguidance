function R = skyBasis(s)
%SKYBASIS Rotation R^c_* of a formation pointing at sky direction s (unit
%   vector in the inertial frame c), before its spin about s: z_* = s,
%   x_* = s x e3 normalized (the e3 x s-free "east" direction; e1 if s is
%   along +-e3), y_* = z_* x x_*. For s = e2 this is Rx(-pi/2), eq. (53).
s = s(:)/norm(s);
x = cross(s, [0; 0; 1]);
if norm(x) < 1e-12, x = [1; 0; 0]; end
x = x/norm(x);
R = [x, cross(s, x), s];
end
