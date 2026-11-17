function c = ocpFunctions(p, ref, withSun)
%OCPFUNCTIONS Mesh, Radau coefficients and the mapped CasADi functions of
%   one spacecraft block (ocpBlock): scaled dynamics (36)-(38) and, if
%   withSun, the Sun cosine (47) at all N*d collocation points. They
%   depend on the physical parameters of p (mass, inertia, SRP), not on
%   its initial state or target point, which enter as NLP parameters.
import casadi.*
c.N = p.ocp.N; c.d = p.ocp.d; c.nx = 13; c.nu = p.nThr;
N = c.N; d = c.d;
c.sigma = ones(1, N)/N;
[c.cRoot, c.C, c.D] = radauCoefficients(d);
c.S = p.ocp.xScale;
c.sNode = [0, cumsum(c.sigma)];
c.sColl = reshape(c.sNode(1:N) + c.cRoot(2:end)'*c.sigma, 1, []);  % column (n-1)*d+j
c.colInt = kron(1:N, ones(1, d));        % interval of each collocation column

% Scaled dynamics dxhat/dt = S^-1 f(t, S xhat, Tmax uhat), eqs. (36)-(38)
pP = pointingParams(p);
t  = SX.sym('t');
xh = SX.sym('x', c.nx);
uh = SX.sym('u', c.nu);
f  = relativeDynamics(t, c.S.*xh, p.Tmax*uh, chebEvalRef(ref, t), pP);
c.F  = Function('F', {t, xh, uh}, {f./c.S});
c.Fm = c.F.map(N*d);
c.Sm = [];
if withSun
    Sig = Function('Sig', {t, xh}, {sunCosine(t, xh(1:4), chebEvalRef(ref, t), pP)});
    c.Sm = Sig.map(N*d);
end
% Physical parameters the functions were built with (setSpacecraft check)
c.physKey = physicsKey(p);
end
