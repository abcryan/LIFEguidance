function s0 = initialState(p)
%INITIALSTATE Truth state matching the paper's rho(0), rhoDot(0).
ic  = p.ic;
Rct = ic.Rte'*ic.Rce;
wct = Rct*ic.wc - ic.wt;
rPcT = Rct*p.rPc;
rct  = ic.rho - rPcT + p.rPt;
rctD = ic.rhoDot - cross(wct, rPcT);
rc = ic.rt + ic.Rte*rct;
vc = ic.vt + ic.Rte*(rctD + cross(ic.wt, rct));
s0 = [ic.Rte(:); ic.wt; ic.rt; ic.vt; ic.Rce(:); ic.wc; rc; vc; zeros(3,1)];
end
