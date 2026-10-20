function w0 = initialGuess(ocp, p)
%INITIALGUESS Sec. 4.3: tau = tf, smooth interpolation of the relative
%   boundary states (cubic Hermite for r, v; slerp for q; linear for w)
%   and T = 0.01 Tmax.
x0 = initialRelativeState(p);
Xall = guessStates([ocp.sNode, ocp.sColl], x0, p.tf);
X  = Xall(:, 1:ocp.N + 1);
XC = Xall(:, ocp.N + 2:end);
U  = p.ocp.uGuess*ones(ocp.nu, ocp.N);
w0 = packW(ocp, 1, X, XC, U);
end

function X = guessStates(s, x0, tf)
r0 = x0(8:10); v0 = x0(11:13)*tf;      % d/ds
h00 = 2*s.^3 - 3*s.^2 + 1; h10 = s.^3 - 2*s.^2 + s;
dh00 = 6*s.^2 - 6*s;       dh10 = 3*s.^2 - 4*s + 1;
r = r0*h00 + v0*h10;
v = (r0*dh00 + v0*dh10)/tf;
q0 = x0(1:4);
ang = acos(min(q0(1), 1));                  % slerp q0 -> [1 0 0 0]
q = (sin((1 - s)*ang).*q0 + sin(s*ang).*[1; 0; 0; 0])/sin(ang);
w = x0(5:7)*(1 - s);
X = [q; w; r; v];
end
