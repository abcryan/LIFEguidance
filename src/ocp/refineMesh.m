function [sigmaNew, marked] = refineMesh(sol, locErr, p)
%REFINEMESH Bisect intervals that limit accuracy or resolution:
%   (a) local discretization error above p.mesh.locErrTol of the gates,
%   (b) thrust switches: intervals whose on/off/saturated pattern differs
%       from a neighbour, and intervals with partially throttled thrusters
%       (an unresolved near-impulsive burn), together with one neighbour on
%       each side.
N = sol.N;
tol = p.mesh.burnTol;
state = zeros(size(sol.T));             % 0 off, 1 partial, 2 saturated
state(sol.T > tol) = 1;
state(sol.T >= p.Tmax*(1 - 1e-3)) = 2;
partial = any(state == 1, 1);
change = [false, any(diff(state, 1, 2) ~= 0, 1)];   % differs from left
switchIv = partial | change | [change(2:end), false];
switchIv = switchIv | [false, switchIv(1:end-1)] | [switchIv(2:end), false];
marked = switchIv | locErr.norm > p.mesh.locErrTol;
if N + nnz(marked) > p.mesh.maxN
    marked(:) = false;                  % cap reached: caller reports it
end
sigmaNew = [];
for n = 1:N
    if marked(n)
        sigmaNew = [sigmaNew, sol.sigma(n)/2, sol.sigma(n)/2]; %#ok<AGROW>
    else
        sigmaNew = [sigmaNew, sol.sigma(n)]; %#ok<AGROW>
    end
end
end
