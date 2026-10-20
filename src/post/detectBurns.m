function burns = detectBurns(sol, prop, p)
%DETECTBURNS Burns = maximal runs of intervals with sum_k T_k > burnTol.
%   Per burn: time span, duration, peak single-thruster thrust, total
%   impulse, and the force and torque integrals in the inertial frame c,
%   int R^c_i f^i dt and int R^c_i tau^i dt (Simpson on the propagated
%   attitude), plus the body-frame integrals as diagnostics.
on = sum(sol.T, 1) > p.mesh.burnTol;
d = diff([false, on, false]);
first = find(d == 1); last = find(d == -1) - 1;
nb = numel(first);
burns = struct('n0', {}, 'n1', {}, 't0', {}, 't1', {}, 'duration', {}, ...
               'peakT', {}, 'impulse', {}, 'Fc', {}, 'Tc', {}, 'Fi', {}, 'Ti', {});
m = prop.nOut;
for b = 1:nb
    n = first(b):last(b);
    B.n0 = first(b); B.n1 = last(b);
    B.t0 = sol.tNode(first(b)); B.t1 = sol.tNode(last(b) + 1);
    B.duration = B.t1 - B.t0;
    B.peakT = max(sol.T(:, n), [], 'all');
    B.impulse = sum(sol.h(n).*sum(sol.T(:, n), 1));
    [f, tau] = thrusterWrench(sol.T(:, n), p);
    B.Fi = f*sol.h(n)'; B.Ti = tau*sol.h(n)';
    B.Fc = zeros(3,1); B.Tc = zeros(3,1);
    for k = n
        idx = (k - 1)*(m - 1) + (1:m);
        tk = prop.t(idx);
        w = simpsonWeights(tk);
        for j = 1:m
            Rci = formationDcm(tk(j), p)*quatToDcm(prop.x(1:4, idx(j)));
            B.Fc = B.Fc + w(j)*Rci*f(:, k - n(1) + 1);
            B.Tc = B.Tc + w(j)*Rci*tau(:, k - n(1) + 1);
        end
    end
    burns(b) = B;
end
end

function w = simpsonWeights(t)
m = numel(t); h = (t(end) - t(1))/(m - 1);
w = h/3*[1, repmat([4 2], 1, (m - 3)/2), 4, 1];
end

function R = formationDcm(t, p)
R = rot1(-pi/2)*rot3(p.phi0 + p.ws*t);
end
