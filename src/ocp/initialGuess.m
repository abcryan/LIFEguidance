function w0 = initialGuess(ocp, p, phiTau, via)
%INITIALGUESS Sec. 5.2: tau = tf, smooth interpolation from x_i,0 to the
%   target state x^tgt(phiTau) of eq. (39) at the guessed arrival phase
%   (cubic Hermite for r, v; slerp for q; linear for w) and T = 0.01 Tmax.
%   Problems 1 and 2: phiTau = phase0 + rate*tf; Problems 3-4: the start's
%   phiTau. via = +1 or -1 (Problem 4, attitude-path class): the attitude
%   is slerped through an intermediate attitude whose z^i points from the
%   anti-Sun direction towards +z_c or -z_c, by p.ocp.p4.viaAntiSunDeg
%   (midway in the field of regard); via = 0 or omitted: direct slerp.
if nargin < 3 || isempty(phiTau) || isempty(ocp.ix.phi)
    phiTau = p.target.phase0 + p.target.rate*p.tf;
end
if nargin < 4, via = 0; end
x0 = initialRelativeState(p);
x1 = targetState(phiTau, p);
s = [ocp.sNode, ocp.sColl];
Xall = guessStates(s, x0, x1, p.tf);
if via ~= 0
    Xall(1:4,:) = viaAttitude(s, x0(1:4), x1(1:4), viaDirection(p, via));
end
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
q = slerp(x0(1:4), x1(1:4), s);
w = x0(5:7)*(1 - s) + x1(5:7)*s;
X = [q; w; r; v];
end

function q = slerp(q0, q1, s)
%SLERP q0 -> nearer of +-q1 at the fractions s (1xn).
if q1'*q0 < 0, q1 = -q1; end
ang = acos(min(q1'*q0, 1));
if ang < 1e-12
    q = repmat(q0, 1, numel(s));
else
    q = (sin((1 - s)*ang).*q0 + sin(s*ang).*q1)/sin(ang);
end
end

function q = viaAttitude(s, q0, q1, dir)
%VIAATTITUDE Slerp q0 -> qMid on s <= 1/2 and qMid -> q1 after, with qMid
%   the direct slerp midpoint rotated minimally so that z^i = dir.
qm = slerp(q0, q1, 0.5);
Rm = quatToDcm(qm);
z = Rm(:,3);
ax = cross(z, dir); sa = norm(ax); ca = z'*dir;
if sa < 1e-12
    Ra = eye(3);
else
    k = ax/sa; ang = atan2(sa, ca);
    Ra = cos(ang)*eye(3) + sin(ang)*skew(k) + (1 - cos(ang))*(k*k');
end
qMid = dcmToQuat(Ra*Rm);
if qMid'*qm < 0, qMid = -qMid; end
q = zeros(4, numel(s));
a = s <= 0.5;
q(:, a) = slerp(q0, qMid, 2*s(a));
q(:, ~a) = slerp(qMid, q1, 2*s(~a) - 1);
end

function dir = viaDirection(p, via)
%VIADIRECTION Unit vector in O*: the anti-Sun direction at tf/2 (from the
%   halo initial position, adequate for a guess) tilted towards via*z_c by
%   p.ocp.p4.viaAntiSunDeg.
fr = frames(p.tf/2, pointingParams(p));  % Rbf = R^b_* (no spin)
a = fr.Rbf'*(p.haloIC(1:3)*p.dse - p.rbs);
a = a/norm(a);
zc = p.target.Rstar'*[0; 0; 1];        % z_c resolved in O*
n = zc - (zc'*a)*a; n = n/norm(n);
dir = cosd(p.ocp.p4.viaAntiSunDeg)*a + via*sind(p.ocp.p4.viaAntiSunDeg)*n;
end
