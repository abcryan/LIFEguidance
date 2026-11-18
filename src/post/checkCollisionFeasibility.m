function chk = checkCollisionFeasibility(p)
%CHECKCOLLISIONFEASIBILITY Before solving Problem 5 (Sec. 5.2, step (i)):
%   the initial formation (initial states, t = 0) and the target formation
%   (independent of phiTau) must satisfy |r_i - r_j| >= dmin for all
%   pairs. Errors if not.
nS = numel(p.fleet.point);
r0 = zeros(3, nS);
for i = 1:nS
    x0 = initialRelativeState(spacecraftParams(p, i));
    r0(:, i) = x0(8:10);
end
P = pyramidPoints(p.target); rT = P(:, p.fleet.point);
pairs = nchoosek(1:nS, 2);
d0 = vecnorm(r0(:, pairs(:,1)) - r0(:, pairs(:,2)));
dT = vecnorm(rT(:, pairs(:,1)) - rT(:, pairs(:,2)));
chk = struct('dmin', p.collision.dmin, 'initialMin', min(d0), 'targetMin', min(dT), ...
             'ok', min(d0) >= p.collision.dmin && min(dT) >= p.collision.dmin);
fprintf(['Collision feasibility: smallest separation %.2f m in the initial formation, ' ...
         '%.2f m in the target formation, dmin %.1f m\n'], chk.initialMin, chk.targetMin, chk.dmin);
if ~chk.ok
    error('LIFE:collisionInfeasible', ['Problem 5 is infeasible: the initial or the target ' ...
          'formation violates the collision constraint (67).']);
end
end
