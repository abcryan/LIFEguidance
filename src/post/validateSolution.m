function res = validateSolution(sol, info, ref, p)
%VALIDATESOLUTION Independent checks of the optimized solution.
%   All reported errors and all animations are based on the forward
%   propagation (propagateRelative), not on the collocation trajectory.
res.sol = sol;
res.info = info;

% Independent forward propagation (dense output for plots and videos)
res.prop = propagateRelative(sol, p, 9);
x = res.prop.x;
res.term = terminalErrors(x(:,end), p);

% Collocation vs propagation at the mesh nodes
res.disc = stateDiscrepancy(res.prop.x(:, res.prop.nodeIdx), sol.X);

% Cross-check with the absolute model (same controls)
tNode = res.prop.t(res.prop.nodeIdx);
res.ab = propagateAbsolute(sol, p, tNode);
res.absDisc = stateDiscrepancy(res.ab.x, res.prop.x(:, res.prop.nodeIdx));
g = p.gate.absRel;
res.absPass = res.absDisc.max.r <= g*p.gate.r && res.absDisc.max.v <= g*p.gate.v ...
    && res.absDisc.max.theta <= g*p.gate.theta && res.absDisc.max.w <= g*p.gate.w;

% Attitude normalization
res.qNormColl = max(abs(vecnorm([sol.X(1:4,:), sol.XC(1:4,:)]) - 1));
res.qNormProp = max(abs(vecnorm(x(1:4,:)) - 1));
Rerr = 0;
for k = 1:size(x, 2)                    % orthogonality of unnormalized R(q)
    q = x(1:4,k); w = q(1); v = q(2:4);
    R = (w^2 - v'*v)*eye(3) + 2*(v*v') + 2*w*skew(v);
    Rerr = max(Rerr, norm(R'*R - eye(3)));
end
res.RorthProp = Rerr;

% Actuation
res.J = sol.J;
res.dm = sol.J/(p.g0*p.Isp);
res.maxThrust = max(sol.T(:));
res.boundViolation = info.boundViolation;
res.thrusterImpulse = sol.T*sol.h';     % 20x1 [N s]
[res.f, res.tau] = thrusterWrench(sol.T, p);
res.cumImpulse = [0, cumsum(sol.h.*sum(sol.T, 1))];
res.burns = detectBurns(sol, res.prop, p);

% Attitude errors along the propagation
[res.thetaAcos, res.thetaAtan] = eigenAngle(x(1:4,:));

% Reference geometry for plots/videos (propagated reference, not the fit)
res.ref = res.prop.ref;
res.Rbf = zeros(3, 3, numel(res.prop.t));
for k = 1:numel(res.prop.t)
    fr = frames(res.prop.t(k), p);
    res.Rbf(:,:,k) = fr.Rbf;
end
haloT = linspace(0, p.haloPeriodDays*86400, 2000);
solFull = ode89(@(~, z) cr3bpNormalized(z, p.eta), [0, p.Om*haloT(end)], p.haloIC, ...
                odeset('RelTol', 1e-12, 'AbsTol', 1e-14));
res.haloFull = deval(solFull, p.Om*haloT);        % normalized
res.refFitErr = ref.fitErr;

res.accepted = info.success && res.term.passAll;
end
