function chk = checkSunFeasibility(p, ref)
%CHECKSUNFEASIBILITY Before solving Problem 2 (Sec. 5.2): the Sun angle of
%   the initial attitude and of the formation attitude (the arrival
%   attitude, z^i = z_*) over [0, tf] must satisfy eq. (48). Errors if not.
p = pointingParams(p);                  % phase-independent (Sec. 4.1)
x0 = initialRelativeState(p);
rbf0 = chebEvalRef(ref, 0);
chk.theta0 = acos(sunCosine(0, x0(1:4), rbf0 + frames(0, p).Rbf*x0(8:10), p));
t = linspace(0, p.tf, 241);
chk.thetaF = acos(sunCosine(t, repmat([1; 0; 0; 0], 1, numel(t)), chebEvalRef(ref, t), p));
inBand = @(th) th >= p.sun.thetaMin & th <= p.sun.thetaMax;
chk.ok = inBand(chk.theta0) && all(inBand(chk.thetaF));
fprintf(['Sun-angle feasibility: anti-Sun angle of n^i at t = 0: %.2f deg, ' ...
         'of the formation over [0, tf]: [%.2f, %.2f] deg, band [%.0f, %.0f] deg\n'], ...
        180 - rad2deg(chk.theta0), 180 - rad2deg(max(chk.thetaF)), 180 - rad2deg(min(chk.thetaF)), ...
        180 - rad2deg(p.sun.thetaMax), 180 - rad2deg(p.sun.thetaMin));
if ~chk.ok
    error('LIFE:sunInfeasible', ['Problem 2 is infeasible: the initial or the ' ...
          'formation attitude violates the Sun-angle constraint (48).']);
end
end
