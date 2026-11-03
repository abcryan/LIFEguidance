function sol = unpackW(ocp, w, p)
%UNPACKW Solution struct in SI units from the decision vector.
w = full(w);
sol.tau   = w(ocp.ix.tau)*p.tf;
sol.sigma = ocp.sigma;
sol.N     = ocp.N;
sol.h     = sol.tau*ocp.sigma;
sol.tNode = sol.tau*ocp.sNode;
sol.tColl = sol.tau*ocp.sColl;
sol.X  = reshape(w(ocp.ix.X), ocp.nx, []).*ocp.S;
sol.XC = reshape(w(ocp.ix.XC), ocp.nx, []).*ocp.S;
sol.U  = reshape(w(ocp.ix.U), ocp.nu, []);
sol.T  = p.Tmax*sol.U;                       % [N], 20xN
sol.J  = sum(sol.h.*sum(sol.T, 1));          % total impulse [N s]
sol.w  = w;
% Arrival phase phiTau and formation phase at t = 0: decision variable in
% Problem 3, prescribed (phi(0) = p.target.phase0) otherwise
if isempty(ocp.ix.phi)
    sol.phase0 = p.target.phase0;
    sol.phiTau = sol.phase0 + p.target.rate*sol.tau;
else
    sol.phiTau = pi*w(ocp.ix.phi);
    sol.phase0 = sol.phiTau - p.target.rate*sol.tau;
end
end
