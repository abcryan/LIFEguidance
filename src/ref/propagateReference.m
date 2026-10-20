function ref = propagateReference(p, tEnd)
%PROPAGATEREFERENCE Unforced formation origin on the L2 halo, eq. (14),
%   integrated in normalized CR3BP units from eq. (51).
%   Returns the ode solution (ref.sol, normalized time), a Chebyshev fit
%   of r^b_bf(t) - r^b_bf(0) on [0, tEnd] for the NLP, and its fit error.
opts = odeset('RelTol', 1e-13, 'AbsTol', 1e-15);
ref.sol = ode89(@(~, z) cr3bpNormalized(z, p.eta), [0, p.Om*tEnd], p.haloIC, opts);
ref.tEnd = tEnd;
ref.eval = @(t) referenceState(ref.sol, t, p);   % SI, t in [s], 6xn

% Chebyshev fit on [0, tEnd] (domain fixed although tau is free)
ref.deg = 14;
nFit = 4*ref.deg;
xi = cos(pi*((nFit:-1:1)' - 0.5)/nFit);          % Chebyshev nodes
tFit = (xi + 1)/2*tEnd;
Z = ref.eval(tFit');
ref.r0 = p.haloIC(1:3)*p.dse;
A = chebBasis(xi, ref.deg);
ref.coef = A\(Z(1:3,:)' - ref.r0');              % (deg+1) x 3

tChk = linspace(0, tEnd, 2001);
Zc = ref.eval(tChk);
fit = chebEvalRef(ref, tChk);
ref.fitErr = max(vecnorm(fit - Zc(1:3,:)));
end

function Z = referenceState(sol, t, p)
z = deval(sol, p.Om*t);
Z = [z(1:3,:)*p.dse; z(4:6,:)*p.dse*p.Om];
end
