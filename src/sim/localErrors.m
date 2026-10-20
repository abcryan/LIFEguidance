function e = localErrors(sol, ref, p)
%LOCALERRORS Per-interval discretization error of the collocation solution:
%   propagate each interval from its collocated start node with the
%   interval's thrust and compare with the collocated end node. Each
%   component is normalized by its arrival gate (59)-(60); e.norm(n) is
%   the largest ratio.
opts = odeset('RelTol', p.prop.RelTol, 'AbsTol', p.prop.AbsTol);
S = p.ocp.xScale;
N = sol.N;
e.r = zeros(1, N); e.v = e.r; e.theta = e.r; e.w = e.r;
for n = 1:N
    t0 = sol.tNode(n); t1 = sol.tNode(n+1);
    zr = deval(ref.sol, p.Om*t0);
    z0 = [sol.X(:,n)./S; zr];
    [~, Z] = ode89(@(tt, zz) relativeRhs(tt, zz, sol.T(:,n), p), [t0, t1], z0, opts);
    d = stateDiscrepancy(S.*Z(end,1:13)', sol.X(:,n+1));
    e.r(n) = d.r; e.v(n) = d.v; e.theta(n) = d.theta; e.w(n) = d.w;
end
e.norm = max([e.r/p.gate.r; e.v/p.gate.v; e.theta/p.gate.theta; e.w/p.gate.w], [], 1);
end
