function zDot = absoluteDynamics(t, z, T, p, c)
%ABSOLUTEDYNAMICS Independent "truth" model: spacecraft and formation
%   origin in the SER frame, eqs. (14)-(15), and inertial spacecraft
%   attitude, eqs. (8)-(9). Positions are shifted by the constant c (the
%   initial formation position) to keep precision.
%   z = [q_ci; w^i_ci; y_i; ydot_i; y_f; ydot_f], r_b = c + y.
q  = z(1:4);  w  = z(5:7);
yi = z(8:10); vi = z(11:13);
yf = z(14:16); vf = z(17:19);
Rci = quatToDcm(q);
Rcb = rot1(p.obliquity)*rot3(p.Om*t);

qDot = 0.5*quatMul(q, [0; w]);          % Rdot^c_i = R^c_i [w^i_ci]x
wDot = p.Jinv*(p.Btau*T - cross(w, p.J*w));

ai = cr3bpAccel(c + yi, vi, p) ...
     + p.beta*(c + yi - p.rbs)/norm(c + yi - p.rbs)^3 ...
     + Rcb'*Rci*(p.Bf*T)/p.m;           % (31)
af = cr3bpAccel(c + yf, vf, p);
zDot = [qDot; wDot; vi; ai; vf; af];
end

function a = cr3bpAccel(r, v, p)
a = 2*p.Om*p.K*v - p.Om^2*(p.K*p.K)*r ...
    - p.mus*(r - p.rbs)/norm(r - p.rbs)^3 - p.mue*(r - p.rbe)/norm(r - p.rbe)^3;
end
