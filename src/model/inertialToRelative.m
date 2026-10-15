function z = inertialToRelative(s, p)
%INERTIALTORELATIVE Map truth state s to the NMPC prediction state z.
Rte = reshape(s(1:9), 3, 3);
wt  = s(10:12);
rt  = s(13:15);
vt  = s(16:18);
Rce = reshape(s(19:27), 3, 3);
wc  = s(28:30);
rc  = s(31:33);
vc  = s(34:36);
xi  = s(37:39);

Rct  = Rte'*Rce;
wct  = Rct*wc - wt;
rPcT = Rct*p.rPc;
rct  = Rte'*(rc - rt);
rctD = Rte'*(vc - vt) - cross(wt, rct);
rho  = rct + rPcT - p.rPt;
rhoD = rctD + cross(wct, rPcT);

z = [Rct(:); wct; rho; rhoD; xi; Rte(:); wt; rt; vt];
end
