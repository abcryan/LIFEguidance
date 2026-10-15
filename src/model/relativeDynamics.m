function zDot = relativeDynamics(z, u, p)
%RELATIVEDYNAMICS NMPC prediction model (15), augmented with target states.
%   z = [R_c/t(:); w^t_c/t; rho; rhoDot; xi;  R_t/e(:); w^t_t/e; r_t; v_t]
%   u = f_ps (24x1), thrust vectors of all thrusters in their F_k.
%   Works for double and CasADi symbolic inputs.

Rct  = reshape(z(1:9), 3, 3);
wct  = z(10:12);
rho  = z(13:15);
rhoD = z(16:18);
Rte  = reshape(z(22:30), 3, 3);
wt   = z(31:33);
rt   = z(34:36);
vt   = z(37:39);

fc   = p.Bf*u;                          % (11)
tauc = p.Btau*u;                        % (12)

% Target: Poisson (1), Euler (2), two-body orbit (6)
RteDot = Rte*skew(wt);
wtDot  = p.JtInv*(skew(p.Jt*wt)*wt);

% Relative attitude: (3) and (4). Since w_c/t is resolved in F_t, Poisson's
% equation is R_c/t_dot = (w^t_c/t)^x R_c/t.
wc     = Rct'*(wct + wt);               % w^c_c/e
wcDot  = p.JcInv*(skew(p.Jc*wc)*wc + tauc);
RctDot = skew(wct)*Rct;
wctDot = Rct*wcDot - wtDot + skew(wct)*(Rct*wc);

% Relative translation of the feature points: (9), (10)
rPcT = Rct*p.rPc;                       % r_pc/c resolved in F_t
rct  = rho - rPcT + p.rPt;              % r_c/t in F_t, from (7)
rctD = rhoD - skew(wct)*rPcT;           % F_t-derivative of r_c/t, from (8)
rc   = rt + Rte*rct;
dAcc = gravity(rc, p.mu) - gravity(rt, p.mu) + Rte*(Rct*fc)/p.mc;
rctDD = Rte'*dAcc - 2*skew(wt)*rctD - skew(wtDot)*rct ...
        - skew(wt)*(skew(wt)*rct);
rhoDD = rctDD + skew(wctDot)*rPcT + skew(wct)*(skew(wct)*rPcT);

zDot = [RctDot(:); wctDot; rhoD; rhoDD; rho;  % (14): xi_dot = rho
        RteDot(:); wtDot; vt; gravity(rt, p.mu)];
end
