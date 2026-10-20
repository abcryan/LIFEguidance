function dz = cr3bpNormalized(z, eta)
%CR3BPNORMALIZED Unforced CR3BP in the SER frame, eq. (14), normalized
%   (length d_se, time 1/Omega). z = [r; v].
r = z(1:3); v = z(4:6);
rs = r - [-eta; 0; 0];
re = r - [1 - eta; 0; 0];
a = [2*v(2) + r(1); -2*v(1) + r(2); 0] ...
    - (1 - eta)*rs/norm(rs)^3 - eta*re/norm(re)^3;
dz = [v; a];
end
