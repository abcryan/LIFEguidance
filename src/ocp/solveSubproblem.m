function r = solveSubproblem(ocp, p, ref, i, phi, tauFix, guess, lam0)
%SOLVESUBPROBLEM One single-spacecraft Problem 3 solve for spacecraft i of
%   the fleet, with the parametric Problem 3 OCP (buildOcp(p, ref, 3,
%   true)).
%   phi: fixed arrival phase [rad], or NaN for free (Problem 3).
%   tauFix: fixed arrival time [s], or NaN for free.
%   guess: initial point w0 (vector), or a struct with fields phi (phase
%   of the target state) and via (attitude class +-1, 0: direct slerp) for
%   the common initial guess (initialGuess).
%   lam0 (optional): initial multipliers (lam_x0, lam_g0) for an OCP built
%   with IPOPT warm-start options.
%   Returns the cost, status, path class (pathClass, from the collocation
%   nodes), the phase sensitivity dJ/dphi = -lam_x(phi)*Jscale/pi
%   [N s/rad] when phi is fixed (envelope theorem; NaN otherwise), and the
%   solution vector with its multipliers.
pk = spacecraftParams(p, i);
o = setSpacecraft(ocp, pk);
if ~isnan(tauFix)
    o.lbx(o.ix.tau) = tauFix/p.tf; o.ubx(o.ix.tau) = tauFix/p.tf;
end
if ~isnan(phi)
    o.lbx(o.ix.phi) = phi/pi; o.ubx(o.ix.phi) = phi/pi;
end
via = NaN;
if isstruct(guess)
    via = guess.via;
    w0 = initialGuess(o, pk, guess.phi, guess.via);
else
    w0 = guess;
end
if nargin < 8, lam0 = []; end
[s, in] = solveOcp(o, w0, pk, [], lam0);
cls = pathClass(s.tNode, s.X(1:4,:), chebEvalRef(ref, s.tNode), pointingParams(pk));
dJ = NaN;
if ~isnan(phi), dJ = -s.lamX(o.ix.phi)*p.ocp.Jscale/pi; end
r = struct('sc', i, 'phiTau', s.phiTau, 'tau', s.tau, 'via', via, 'J', s.J, ...
    'status', string(in.status), 'success', in.success, 'iter', in.iter, ...
    'time', in.wallTime, 'cls', cls, 'dJdphi', dJ, 'w', s.w, 'lamX', s.lamX, ...
    'lamG', s.lamG);
end
