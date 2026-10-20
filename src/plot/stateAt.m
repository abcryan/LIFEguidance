function [x, T, k] = stateAt(res, t)
%STATEAT Propagated relative state and thrust at physical times t (1xn):
%   linear interpolation of the dense propagation output (quaternion
%   renormalized), piecewise-constant thrust of the optimal controls.
tp = res.prop.t;
t = min(max(t, tp(1)), tp(end));       % sum(sigma) = 1 only to round-off
x = interp1(tp, res.prop.x', t, 'linear')';
x(1:4,:) = x(1:4,:)./vecnorm(x(1:4,:));
k = discretize(t, res.sol.tNode);
k(isnan(k)) = res.sol.N;
T = res.sol.T(:, k);
end
