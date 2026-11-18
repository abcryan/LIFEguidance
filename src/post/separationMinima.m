function S = separationMinima(res, p, pairs)
%SEPARATIONMINIMA Every local minimum of the distance d_ij(t) between the
%   forward-propagated spacecraft (Problem 5 verification (iii), Sec. 5.2).
%   res{i}: validated spacecraft (validateFleet) on the common propagation
%   grid (shared tau and mesh). Candidates: sign changes of
%   d/dt d_ij^2 = 2 dr'*dv from - to + between output samples, and the two
%   ends. Each interior candidate is located exactly: both spacecraft are
%   integrated again over that control interval from their propagated
%   state at its start node (ode89, same tolerances and thrust as
%   propagateRelative), and the root is found on the integrator's dense
%   output (deval, closestApproach).
%   S(k): pair, t and d of every local minimum, and the global minimum.
t = res{1}.prop.t;
for i = 2:numel(res)
    assert(isequal(res{i}.prop.t, t), 'Propagation grids differ.');
end
nOut = res{1}.prop.nOut; nodeIdx = res{1}.prop.nodeIdx;
Sx = p.ocp.xScale;
opts = odeset('RelTol', p.prop.RelTol, 'AbsTol', p.prop.AbsTol);
cache = containers.Map('KeyType', 'char', 'ValueType', 'any');
S = struct('pair', {}, 't', {}, 'd', {}, 'tMin', {}, 'dMin', {});
for k = 1:size(pairs, 1)
    i = pairs(k, 1); j = pairs(k, 2);
    dX = res{i}.prop.xP(8:13, :) - res{j}.prop.xP(8:13, :);
    gs = sum(dX(1:3, :).*dX(4:6, :), 1);
    dd = vecnorm(dX(1:3, :));
    cand = find(gs(1:end-1) < 0 & gs(2:end) >= 0);
    tm = []; dm = [];
    for c = cand
        n = floor((c - 1)/(nOut - 1)) + 1;          % control interval of sample c
        si = intervalSol(res{i}, n, nodeIdx, Sx, p, opts, cache, sprintf('%d-%d', i, n));
        sj = intervalSol(res{j}, n, nodeIdx, Sx, p, opts, cache, sprintf('%d-%d', j, n));
        rv = @(tt) Sx(8:13).*(deval(si, tt, 8:13) - deval(sj, tt, 8:13));
        [tc, dc] = closestApproach(rv, t(c), t(c + 1));
        tm(end + 1) = tc; dm(end + 1) = dc; %#ok<AGROW>
    end
    if gs(1) >= 0, tm = [t(1), tm]; dm = [dd(1), dm]; end          %#ok<AGROW>
    if gs(end) <= 0, tm = [tm, t(end)]; dm = [dm, dd(end)]; end    %#ok<AGROW>
    [dMin, km] = min(dm);
    S(k) = struct('pair', pairs(k, :), 't', tm, 'd', dm, 'tMin', tm(km), 'dMin', dMin);
end
end

function sol = intervalSol(r, n, nodeIdx, Sx, p, opts, cache, key)
% ode89 solution of spacecraft r over control interval n from its
% propagated state at the start node (as in propagateRelative).
if isKey(cache, key), sol = cache(key); return, end
kk = nodeIdx(n);
pP = pointingParams(r.p);
z0 = [r.prop.xP(:, kk)./Sx; r.prop.ref(1:3, kk)/p.dse; r.prop.ref(4:6, kk)/(p.dse*p.Om)];
T = r.sol.T(:, n);
sol = ode89(@(tt, zz) relativeRhs(tt, zz, T, pP), [r.sol.tNode(n), r.sol.tNode(n+1)], z0, opts);
cache(key) = sol; %#ok<NASGU> (handle: stored in the caller's map)
end
