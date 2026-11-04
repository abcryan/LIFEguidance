function x = targetState(phiTau, p)
%TARGETSTATE Target state x_i^tgt(phiTau) in the pointing frame O*,
%   eq. (39): R*_i = Rz(phiTau), w*_*i = phidot e3, r*_*i = Rz(phiTau) r_fp,
%   rdot*_*i = phidot e3 x Rz(phiTau) r_fp. Double or CasADi phiTau.
Rz = rot3(phiTau);
w = [0; 0; p.target.rate];
r = Rz*p.rfp;
x = [cos(phiTau/2); 0; 0; sin(phiTau/2); w; r; cross(w, r)];
end
