function x0 = initialRelativeState(p)
%INITIALRELATIVESTATE x_i(0) from the initial virtual formation p.init.
%   The initial formation exists only at t = 0. It shares origin (position
%   and velocity) with the target formation, points at p.init.sky and spins
%   at p.init.rate about its own z axis. The spacecraft sits at its point
%   p.init.point, with body axes aligned with the initial frame, and
%   co-rotates with it, so
%     R^c_i = R^c_f0,  w^i_ci = rate*e3,  d/dt|_c r_fi = w_cf0 x r_fi.
spec = p.init;
Rcf0 = formationAttitude(spec, 0);
Rcf  = formationAttitude(p.target, 0);
Rf_f0 = Rcf'*Rcf0;                         % initial frame -> target frame
Rfi = Rf_f0;                               % body axes = initial-frame axes
wfi = Rfi*[0; 0; spec.rate] - [0; 0; p.target.rate];      % (11)
P = pyramidPoints(spec);
rfi = Rf_f0*P(:, spec.point);              % r^f_fi
wRel = Rf_f0*[0; 0; spec.rate] - [0; 0; p.target.rate];   % w^f_cf0 - w^f_cf
vfi = cross(wRel, rfi);                    % f-frame derivative of r_fi
x0 = [dcmToQuat(Rfi); wfi; rfi - p.rfp; vfi];
end
