function fr = frames(t, p)
%FRAMES Prescribed SER and formation attitudes at time t [s].
%   Eq. (52), the target formation of p.target (eq. (53) by default) and
%   (23)-(27). Fields:
%     Rcb, Rcf  R^c_b, R^c_f (b, f to the inertial ICR frame c)
%     Rbf       R^b_f = (R^c_b)' R^c_f
%     wcf, wcb  w^f_cf, w^b_cb
%     wbf       w^f_bf, eq. (24)
%     wbfDot    d/dt w^f_bf, eq. (27) (uniform SER rotation, wcf constant)
%   Double or CasADi t.
fr.Rcb = rot1(p.obliquity)*rot3(p.Om*t);
fr.Rcf = formationAttitude(p.target, t);
fr.Rbf = fr.Rcb'*fr.Rcf;
fr.wcb = [0; 0; p.Om];
fr.wcf = [0; 0; p.ws];
wcbF = fr.Rbf'*fr.wcb;                  % w_cb resolved in f
fr.wbf = fr.wcf - wcbF;
fr.wbfDot = skew(fr.wbf)*wcbF;
end
