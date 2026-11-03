function x = shiftPhase(x0, dphi, p)
%SHIFTPHASE Relative state of the same physical spacecraft state w.r.t. the
%   target formation advanced in phase by dphi about z_f, eq. (47):
%   R^c_f -> R^c_f Rz(dphi), so every quantity resolved in O^f is rotated
%   by Rz(-dphi) and R^f_i -> Rz(-dphi) R^f_i. w_cf = rate*e3 is unchanged.
%   x0 is the relative state (13x1) at the reference phase, e.g. from
%   initialRelativeState(p). Double or CasADi dphi (closed form, unlike
%   initialRelativeState with a shifted p.target.phase0, which it equals).
R = rot3(-dphi);
q = quatMul([cos(dphi/2); 0; 0; -sin(dphi/2)], x0(1:4));
x = [q; R*x0(5:7); R*(x0(8:10) + p.rfp) - p.rfp; R*x0(11:13)];
end
