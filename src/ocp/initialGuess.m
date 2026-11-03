function w0 = initialGuess(ocp, p, phiTau)
%INITIALGUESS Sec. 5.2: tau = tf, smooth interpolation of the relative
%   boundary states (cubic Hermite for r, v; slerp for q; linear for w)
%   and T = 0.01 Tmax. Problem 3: the boundary states are those of the
%   start's arrival phase phiTau, x_i,0(phiTau, tf) of eq. (47).
pS = p;
if nargin < 3, phiTau = []; end
if ~isempty(ocp.ix.phi)
    pS.target.phase0 = phiTau - p.target.rate*p.tf;
end
x0 = initialRelativeState(pS);
Xall = guessStates([ocp.sNode, ocp.sColl], x0, p.tf);
X  = Xall(:, 1:ocp.N + 1);
XC = Xall(:, ocp.N + 2:end);
U  = p.ocp.uGuess*ones(ocp.nu, ocp.N);
w0 = packW(ocp, 1, X, XC, U, phiTau);
end

function X = guessStates(s, x0, tf)
r0 = x0(8:10); v0 = x0(11:13)*tf;      % d/ds
h00 = 2*s.^3 - 3*s.^2 + 1; h10 = s.^3 - 2*s.^2 + s;
dh00 = 6*s.^2 - 6*s;       dh10 = 3*s.^2 - 4*s + 1;
r = r0*h00 + v0*h10;
v = (r0*dh00 + v0*dh10)/tf;
q0 = x0(1:4);
q1 = [sign(q0(1)) + (q0(1) == 0); 0; 0; 0]; % nearer of +-identity (q_v = 0)
ang = acos(min(q1'*q0, 1));                 % slerp q0 -> q1
q = (sin((1 - s)*ang).*q0 + sin(s*ang).*q1)/sin(ang);
w = x0(5:7)*(1 - s);
X = [q; w; r; v];
end
