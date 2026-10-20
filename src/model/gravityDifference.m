function da = gravityDifference(R, d, mu)
%GRAVITYDIFFERENCE g(R + d) - g(R) for g(r) = -mu*r/|r|^3, without the
%   cancellation of a direct subtraction (Battin's f(q) form).
%   R: position of the reference point relative to the attractor, d: offset.
%   Double or CasADi.
q  = (d'*(d + 2*R))/(R'*R);
fq = q*(3 + 3*q + q^2)/(1 + (1 + q)^1.5);   % (1+q)^(3/2) - 1
rho = R + d;
da = -mu*(d - fq*R)/(rho'*rho)^1.5;
end
