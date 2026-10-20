function z = relToAbs(t, x, yf, vf, p)
%RELTOABS Absolute state of absoluteDynamics from the relative state (33)
%   and the (shifted) formation origin state yf, vf in SER.
fr = frames(t, p);
qci = quatMul(formationQuat(t, p), x(1:4));
w = quatToDcm(x(1:4))'*(x(5:7) + fr.wcf);
rfi = x(8:10) + p.rfp;
yi = yf + fr.Rbf*rfi;
vi = vf + fr.Rbf*(x(11:13) + cross(fr.wbf, rfi));
z = [qci; w; yi; vi; yf; vf];
end
