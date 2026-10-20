function w0 = warmStart(ocp, prev, p)
%WARMSTART Interpolate a previous solution (any mesh) onto a new mesh in s.
sPrev = [prev.tNode, prev.tColl]/prev.tau;
[sPrev, k] = unique(sPrev);
Xprev = [prev.X, prev.XC];
Xprev = Xprev(:, k);
sAll = [ocp.sNode, ocp.sColl];
Xall = interp1(sPrev, Xprev', sAll, 'pchip')';
X  = Xall(:, 1:ocp.N + 1);
XC = Xall(:, ocp.N + 2:end);
% Piecewise-constant controls, sampled at the new interval midpoints
sMid = (ocp.sNode(1:end-1) + ocp.sNode(2:end))/2;
edges = prev.tNode/prev.tau;
n = discretize(sMid, edges);
U = prev.U(:, n);
w0 = packW(ocp, prev.tau/p.tf, X, XC, U);
end
