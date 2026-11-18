function [tc, dc] = closestApproach(rv, ta, tb)
%CLOSESTAPPROACH Time and distance of the local minimum of |dr(t)| in
%   [ta, tb], where rv(t) = [dr; dv] (6x1) is the relative position and
%   velocity of two spacecraft: the root of d/dt |dr|^2 = 2 dr'*dv, which
%   must change sign from - to + in [ta, tb]; otherwise the smaller end
%   value is returned (a minimum at an end of the bracket).
g = @(t) dotRv(rv(t));
ga = g(ta); gb = g(tb);
if ga < 0 && gb > 0
    tc = fzero(g, [ta tb], optimset('TolX', 1e-9*max(1, abs(tb))));
else
    da = norm(head3(rv(ta))); db = norm(head3(rv(tb)));
    if da <= db, tc = ta; else, tc = tb; end
end
dc = norm(head3(rv(tc)));
end

function v = dotRv(x)
v = x(1:3)'*x(4:6);
end

function r = head3(x)
r = x(1:3);
end
