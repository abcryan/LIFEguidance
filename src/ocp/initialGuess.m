function w0 = initialGuess(ocp, p, phiTau)
%INITIALGUESS Sec. 5.2: tau = tf, smooth interpolation from x_i,0 to the
%   target state x^tgt(phiTau) of eq. (39) at the guessed arrival phase
%   (cubic Hermite for r, v; slerp for q; linear for w) and T = 0.01 Tmax.
%   Problems 1 and 2: phiTau = phase0 + rate*tf; Problem 3: the start's
%   phiTau.
if nargin < 3 || isempty(ocp.ix.phi)
    phiTau = p.target.phase0 + p.target.rate*p.tf;
end
x0 = initialRelativeState(p);
x1 = targetState(phiTau, p);
Xall = guessStates([ocp.sNode, ocp.sColl], x0, x1, p.tf);
X  = Xall(:, 1:ocp.N + 1);
XC = Xall(:, ocp.N + 2:end);
U  = p.ocp.uGuess*ones(ocp.nu, ocp.N);
if isempty(ocp.ix.phi), phiTau = []; end
w0 = packW(ocp, 1, X, XC, U, phiTau);
end

function X = guessStates(s, x0, x1, tf)
r0 = x0(8:10); v0 = x0(11:13)*tf;      % d/ds
r1 = x1(8:10); v1 = x1(11:13)*tf;
h00 = 2*s.^3 - 3*s.^2 + 1; h10 = s.^3 - 2*s.^2 + s;
h01 = -2*s.^3 + 3*s.^2;    h11 = s.^3 - s.^2;
dh00 = 6*s.^2 - 6*s;       dh10 = 3*s.^2 - 4*s + 1;
dh01 = -6*s.^2 + 6*s;      dh11 = 3*s.^2 - 2*s;
r = r0*h00 + v0*h10 + r1*h01 + v1*h11;
v = (r0*dh00 + v0*dh10 + r1*dh01 + v1*dh11)/tf;
q0 = x0(1:4);
q1 = x1(1:4);
if q1'*q0 < 0, q1 = -q1; end            % nearer of +-q^tgt
ang = acos(min(q1'*q0, 1));             % slerp q0 -> q1
if ang < 1e-12
    q = repmat(q0, 1, numel(s));
else
    q = (sin((1 - s)*ang).*q0 + sin(s*ang).*q1)/sin(ang);
end
w = x0(5:7)*(1 - s) + x1(5:7)*s;
X = [q; w; r; v];
end
