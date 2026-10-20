function ab = propagateAbsolute(sol, p, tOut)
%PROPAGATEABSOLUTE Cross-check: the same controls applied to the absolute
%   model (spacecraft and formation origin in SER, inertial attitude),
%   mapped back to the relative state (33) at times tOut.
%   Positions are shifted by c = r_bf(0) to retain precision.
c = p.haloIC(1:3)*p.dse;
vf0 = p.haloIC(4:6)*p.dse*p.Om;
z = relToAbs(0, initialRelativeState(p), zeros(3,1), vf0, p);
opts = odeset('RelTol', 1e-13, ...
              'AbsTol', [1e-15*ones(4,1); 1e-16*ones(3,1); ...
                         1e-7*ones(3,1); 1e-11*ones(3,1); ...
                         1e-7*ones(3,1); 1e-11*ones(3,1)]);
ab.t = tOut;
ab.x = zeros(13, numel(tOut));
ab.z = zeros(19, numel(tOut));
for n = 1:sol.N
    t0 = sol.tNode(n); t1 = sol.tNode(n+1);
    inside = tOut > t0 & tOut <= t1;
    if n == 1, inside = inside | tOut == t0; end
    tn = unique([t0, tOut(inside), t1]);
    if numel(tn) == 2, tn = [t0, (t0 + t1)/2, t1]; end
    [tt, Zn] = ode89(@(tt, zz) absoluteDynamics(tt, zz, sol.T(:,n), p, c), tn, z, opts);
    for k = find(inside)
        j = find(tt == tOut(k), 1);
        ab.z(:,k) = Zn(j,:)';
        ab.x(:,k) = absToRel(tOut(k), Zn(j,:)', p);
    end
    z = Zn(end,:)';
end
end
