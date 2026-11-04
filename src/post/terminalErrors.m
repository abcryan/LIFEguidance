function e = terminalErrors(xP, phiTau, p)
%TERMINALERRORS Arrival checks of Sec. 5.2 on the pointing-frame state xP
%   (13x1) at arrival: deviations of r, rdot and w from x^tgt(phiTau),
%   eq. (39), and the eigenaxis angle xi of R^f_i = Rz(phiTau)' R*_i.
xt = targetState(phiTau, p);
qfi = quatMul([cos(phiTau/2); 0; 0; -sin(phiTau/2)], xP(1:4));
e.r = norm(xP(8:10) - xt(8:10));
e.v = norm(xP(11:13) - xt(11:13));
[e.thetaAcos, e.thetaAtan] = eigenAngle(qfi);
e.w = norm(xP(5:7) - xt(5:7));
e.qNormErr = abs(norm(xP(1:4)) - 1);
e.pass.r = e.r <= p.gate.r;
e.pass.v = e.v <= p.gate.v;
e.pass.theta = e.thetaAcos <= p.gate.theta;
e.pass.w = e.w <= p.gate.w;
e.passAll = e.pass.r && e.pass.v && e.pass.theta && e.pass.w;
end
