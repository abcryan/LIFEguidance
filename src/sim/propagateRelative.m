function prop = propagateRelative(sol, p, nOut)
%PROPAGATERELATIVE Independent forward propagation of the optimized
%   piecewise-constant thrust: pointing-frame model (36)-(38)
%   co-integrated with the reference ODE (14) (not the Chebyshev fit used
%   in the NLP). prop.xP is the pointing-frame state, prop.x the
%   target-relative formation-frame state (toFormation, phase of p).
%   ode89 restarted at every control interval, RelTol/AbsTol of Sec. 4.3
%   on states scaled as in the NLP (reference in CR3BP units).
%   nOut: output samples per interval (>= 2, includes both ends).
if nargin < 3, nOut = 9; end
S = p.ocp.xScale;
opts = odeset('RelTol', p.prop.RelTol, 'AbsTol', p.prop.AbsTol);
pP = pointingParams(p);
z = [initialRelativeState(p)./S; p.haloIC];
N = sol.N;
t = zeros(1, N*(nOut - 1) + 1);
Z = zeros(19, numel(t));
T = zeros(p.nThr, numel(t));
Z(:,1) = z; T(:,1) = sol.T(:,1);
for n = 1:N
    tn = linspace(sol.tNode(n), sol.tNode(n+1), nOut);
    [~, Zn] = ode89(@(tt, zz) relativeRhs(tt, zz, sol.T(:,n), pP), tn, z, opts);
    if nOut == 2, Zn = Zn([1 end],:); end   % two-point tspan returns all steps
    k = (n - 1)*(nOut - 1) + (2:nOut);
    t(k) = tn(2:end);
    Z(:,k) = Zn(2:end,:)';
    T(:,k) = repmat(sol.T(:,n), 1, nOut - 1);
    z = Zn(end,:)';
end
prop.t = t;
prop.xP = Z(1:13,:).*S;                 % SI pointing-frame state
prop.x = toFormation(t, prop.xP, p);    % target-relative, formation frame
prop.ref = [Z(14:16,:)*p.dse; Z(17:19,:)*p.dse*p.Om];   % r_bf, v_bf [m, m/s]
prop.T = T;                             % thrust held on [t(k-1), t(k)]
prop.nOut = nOut;
prop.nodeIdx = 1:(nOut - 1):numel(t);
end

