function sun = sunAngle(res, ref, p)
%SUNANGLE Sun angle theta_i of n^i, eqs. (47)-(48), for any problem.
%   Along the propagation (res.prop) with the exact Sun direction from the
%   spacecraft: theta_i, the anti-Sun angle 180 - theta_i, the bound
%   violation and the intervals outside [thetaMin, thetaMax] by more than
%   the gate p.gate.sunDeg (smaller inter-sample dips of a constrained
%   solution show in the max violation only), the error of
%   the formation-centre direction used in the NLP, and the azimuth of n^i
%   about the Sun line (sky plot). At the collocation points, the NLP form
%   (formation centre, fitted reference). For Problem 2 also the points
%   where the constraint is active.
t = res.prop.t; x = res.prop.x; n = numel(t);
rbf = res.prop.ref(1:3,:);
rbi = zeros(3, n); az = zeros(1, n);
for k = 1:n
    rbi(:,k) = rbf(:,k) + res.Rbf(:,:,k)*(x(8:10,k) + p.rfp);
    a = rbi(:,k) - p.rbs; a = a/norm(a);                 % anti-Sun direction
    u1 = [0; 0; 1] - a(3)*a; u1 = u1/norm(u1);           % ecliptic north
    nb = res.Rbf(:,:,k)*quatToDcm(x(1:4,k))*p.sun.n;
    az(k) = atan2(nb'*cross(a, u1), nb'*u1);
end
sig = sunCosine(t, x(1:4,:), rbi, p);
sigC = sunCosine(t, x(1:4,:), rbf, p);
sun.t = t;
sun.theta = acos(min(max(sig, -1), 1));
sun.antiDeg = 180 - rad2deg(sun.theta);
sun.azimuth = az;
sun.centreErr = max(abs(sun.theta - acos(min(max(sigC, -1), 1))));
[sun.viol, sun.maxViolDeg] = violationDeg(sun.theta, p);
out = sun.viol > p.gate.sunDeg;
sun.timeOutside = sum(diff(t).*(out(1:end-1) | out(2:end)));
sun.intervals = runs(t, out, sun.antiDeg);
sun.pass = sun.maxViolDeg <= p.gate.sunDeg;

% NLP form at the collocation points
sol = res.sol;
sigN = sunCosine(sol.tColl, sol.XC(1:4,:), chebEvalRef(ref, sol.tColl), pointingParams(p));
[~, sun.collMaxViolDeg] = violationDeg(acos(min(max(sigN, -1), 1)), p);

% Active constraint (Problem 2): collocation points within the gate of a
% bound. With the attitude a path constraint of high relative degree, these
% are isolated touch points rather than boundary arcs. The CasADi
% multiplier there is negative at the lower bound of sigma (thetaMax,
% anti-Sun 50 deg) and positive at the upper one (thetaMin, anti-Sun 85).
sun.active = struct('t', {}, 'antiDeg', {}, 'bound', {}, 'lambda', {});
sun.activeBound = "";
if isfield(sol, 'lamSun') && ~isempty(sol.lamSun)
    thN = rad2deg(acos(min(max(sigN, -1), 1)));
    lo = abs(thN - rad2deg(p.sun.thetaMax)) < p.gate.sunDeg;
    hi = abs(thN - rad2deg(p.sun.thetaMin)) < p.gate.sunDeg;
    for k = find(lo | hi)
        sun.active(end + 1) = struct('t', sol.tColl(k), 'antiDeg', 180 - thN(k), ...
            'bound', 180 - rad2deg(p.sun.thetaMax*lo(k) + p.sun.thetaMin*hi(k)), ...
            'lambda', sol.lamSun(k));
    end
    sides = ["anti-Sun 50 deg", "anti-Sun 85 deg"];
    sun.activeBound = strjoin(sides([any(lo), any(hi)]), ", ");
end
end

function [v, vmax] = violationDeg(theta, p)
v = rad2deg(max(max(p.sun.thetaMin - theta, theta - p.sun.thetaMax), 0));
vmax = max([v, 0]);
end

function iv = runs(t, mask, antiDeg)
%RUNS Maximal runs of mask over the samples t: start, end and the extreme
%   anti-Sun angle in the run.
iv = struct('t0', {}, 't1', {}, 'antiDegMin', {}, 'antiDegMax', {});
if isempty(mask), return; end
d = diff([false, mask(:)', false]);
first = find(d == 1); last = find(d == -1) - 1;
for b = 1:numel(first)
    k = first(b):last(b);
    iv(b).t0 = t(first(b));
    iv(b).t1 = t(min(last(b) + 1, numel(t)));
    iv(b).antiDegMin = min(antiDeg(k));
    iv(b).antiDegMax = max(antiDeg(k));
end
end
