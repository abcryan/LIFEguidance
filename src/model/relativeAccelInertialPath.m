function rDD = relativeAccelInertialPath(t, x, T, rbf, p)
%RELATIVEACCELINERTIALPATH Check of (19)-(22) by an independent route:
%   difference the inertial accelerations of i and f, then transform to
%   the formation frame, which spins at the constant rate w^f_cf:
%   rddot^f = R^f_b G^b - 2 w_cf x rdot^f - w_cf x (w_cf x r^f_fi).
q = x(1:4); r = x(8:10); v = x(11:13);
fr = frames(t, p);
Rfi = quatToDcm(q);
rfi = r + p.rfp;
rbi = rbf + fr.Rbf*rfi;
g = @(r, c, mu) -mu*(r - c)/norm(r - c)^3;
G = g(rbi, p.rbs, p.mus) - g(rbf, p.rbs, p.mus) ...
    + g(rbi, p.rbe, p.mue) - g(rbf, p.rbe, p.mue) ...
    + p.beta*(rbi - p.rbs)/norm(rbi - p.rbs)^3 + fr.Rbf*Rfi*(p.Bf*T)/p.m;
rDD = fr.Rbf'*G - 2*cross(fr.wcf, v) - cross(fr.wcf, cross(fr.wcf, rfi));
end
