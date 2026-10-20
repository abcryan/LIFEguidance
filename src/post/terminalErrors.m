function e = terminalErrors(x, p)
%TERMINALERRORS Arrival checks (59)-(60) on a relative state x (13x1).
e.r = norm(x(8:10));
e.v = norm(x(11:13));
[e.thetaAcos, e.thetaAtan] = eigenAngle(x(1:4));
e.w = norm(x(5:7));
e.qNormErr = abs(norm(x(1:4)) - 1);
e.pass.r = e.r <= p.gate.r;
e.pass.v = e.v <= p.gate.v;
e.pass.theta = e.thetaAcos <= p.gate.theta;
e.pass.w = e.w <= p.gate.w;
e.passAll = e.pass.r && e.pass.v && e.pass.theta && e.pass.w;
end
