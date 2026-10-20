function [xDot, out] = relativeDynamics(t, x, T, rbf, p)
%RELATIVEDYNAMICS Spacecraft state relative to the target point p on the
%   rotating virtual formation, eqs. (10), (13), (19)-(22), (29)-(31).
%   x   = [q; w^f_fi; r^f_pi; rdot^f_pi] (13x1, SI), q ~ R^f_i
%         (Hamilton, scalar first)
%   T   = thrust magnitudes (20x1) [N]
%   rbf = formation origin in the SER frame [m] (reference trajectory)
%   Works for double and CasADi inputs.
%   out (optional): body-frame force/torque and frames for diagnostics.
q = x(1:4); w = x(5:7); r = x(8:10); v = x(11:13);
fr = frames(t, p);
Rfi = quatToDcm(q);

% Control wrench in O^i, (29)-(30)
fi = p.Bf*T;
taui = p.Btau*T;

% Relative attitude: quaternion form of (10), and (13)
qDot = 0.5*[-w'*q(2:4); q(1)*w + skew(w)*q(2:4)];   % [0; w] (x) q
wciF = w + fr.wcf;                      % w^f_ci, from (11)
wci  = Rfi'*wciF;                       % w^i_ci
wciDot = p.Jinv*(taui - skew(wci)*(p.J*wci));       % Euler (9)
wDot = Rfi*wciDot + skew(w)*(Rfi*wci);  % (13) with wdot^f_cf = 0

% Relative translation, (19)-(22). d = r_bi - r_bf.
rfi = r + p.rfp;                        % r^f_fi, from (16)
dB  = fr.Rbf*rfi;                       % resolved in b
dvB = fr.Rbf*(v + skew(fr.wbf)*rfi);    % rdot^b_bi - rdot^b_bf
rbi = rbf + dB;
srp = p.beta*(rbi - p.rbs)/((rbi - p.rbs)'*(rbi - p.rbs))^1.5;
ddB = 2*p.Om*p.K*dvB - p.Om^2*(p.K*p.K)*dB ...       % (15) - (14)
      + gravityDifference(rbf - p.rbs, dB, p.mus) ...
      + gravityDifference(rbf - p.rbe, dB, p.mue) ...
      + srp + fr.Rbf*(Rfi*fi)/p.m;      % a_ctrl via (31)
Rfb = fr.Rbf';
rDD = Rfb*ddB - 2*skew(fr.wbf)*(Rfb*dvB) - skew(fr.wbfDot)*rfi ...
      + skew(fr.wbf)*(skew(fr.wbf)*rfi);

xDot = [qDot; wDot; v; rDD];
if nargout > 1
    out.fi = fi; out.taui = taui; out.Rfi = Rfi; out.fr = fr; out.srp = srp;
end
end
