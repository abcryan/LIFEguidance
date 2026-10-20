function q = formationQuat(t, p)
%FORMATIONQUAT Quaternion of R^c_f(t) = Rx(-pi/2) Rz(phi0 + ws*t), eq. (53),
%   continuous in t.
a = -pi/2; b = p.phi0 + p.ws*t;
q = quatMul([cos(a/2); sin(a/2); 0; 0], [cos(b/2); 0; 0; sin(b/2)]);
end
