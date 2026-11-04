function x0 = initialRelativeState(p)
%INITIALRELATIVESTATE x_i,0 in the pointing frame O* of the target
%   formation, from the initial virtual formation p.init. The initial
%   formation exists only at t = 0, shares the target origin (position and
%   velocity), points at p.init.sky and spins at p.init.rate about its own
%   z axis. The spacecraft sits at its point p.init.point, with body axes
%   aligned with the initial frame, and co-rotates with it, so
%     R^c_i = R^c_f0,  w^i_ci = rate*e3,  d/dt|_c r_fi = w_cf0 x r_fi.
%   Independent of the target phase (Sec. 4.1, boundary conditions).
spec = p.init;
R = p.target.Rstar'*formationAttitude(spec, 0);   % R*_f0 = R*_i
w = R*[0; 0; spec.rate];                % w*_*i = w*_ci (w_c* = 0)
P = pyramidPoints(spec);
r = R*P(:, spec.point);                 % r*_*i
x0 = [dcmToQuat(R); w; r; cross(w, r)];
end
