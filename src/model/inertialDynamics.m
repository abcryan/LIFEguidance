function sDot = inertialDynamics(s, u, p)
%INERTIALDYNAMICS Truth model: target and chaser rigid bodies in ECI.
%   s = [R_t/e(:); w^t_t/e; r_t; v_t;  R_c/e(:); w^c_c/e; r_c; v_c;  xi]
%   Eqs. (1), (2), (6) with zero perturbations, plus the integrator (14).

Rte = reshape(s(1:9), 3, 3);
wt  = s(10:12);
rt  = s(13:15);
vt  = s(16:18);
Rce = reshape(s(19:27), 3, 3);
wc  = s(28:30);
rc  = s(31:33);
vc  = s(34:36);

fc   = p.Bf*u;
tauc = p.Btau*u;

RteDot = Rte*skew(wt);
wtDot  = p.JtInv*(skew(p.Jt*wt)*wt);
RceDot = Rce*skew(wc);
wcDot  = p.JcInv*(skew(p.Jc*wc)*wc + tauc);

rho = Rte'*(rc - rt) + Rte'*Rce*p.rPc - p.rPt;

% Integral action (14), with conditional integration as anti-windup
xiDot = rho*(norm(rho) <= p.mpc.xiRadius);

sDot = [RteDot(:); wtDot; vt; gravity(rt, p.mu);
        RceDot(:); wcDot; vc; gravity(rc, p.mu) + Rce*fc/p.mc;
        xiDot];
end
