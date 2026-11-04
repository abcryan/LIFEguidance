function xf = toFormation(t, xP, p)
%TOFORMATION Target-relative state in the spinning formation frame O^f,
%   [q_fi; w^f_fi; r^f_pi; rdot^f_pi] (the layout of the formation-frame
%   relative model), from the pointing-frame state xP (13xn) at times t
%   (1xn), with R*_f = Rz(phi(t)) and phi(t) = p.target.phase0 + rate*t,
%   eq. (33).
rate = p.target.rate;
w = [0; 0; rate];
n = size(xP, 2);
xf = zeros(13, n);
for k = 1:n
    phi = p.target.phase0 + rate*t(k);
    Rt = rot3(phi)';                    % R^f_*
    r = xP(8:10,k);
    xf(:,k) = [quatMul([cos(phi/2); 0; 0; -sin(phi/2)], xP(1:4,k));
               Rt*(xP(5:7,k) - w);
               Rt*r - p.rfp;
               Rt*(xP(11:13,k) - cross(w, r))];
end
end
