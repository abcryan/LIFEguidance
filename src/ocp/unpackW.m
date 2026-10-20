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
end
