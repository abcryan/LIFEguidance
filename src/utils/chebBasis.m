function T = chebBasis(x, deg)
%CHEBBASIS Chebyshev polynomials T_0..T_deg at x in [-1, 1]. Double or CasADi.
%   Double input: x is a column vector, T is numel(x) x (deg+1).
%   CasADi input: x is scalar, T is 1 x (deg+1).
if isa(x, 'double')
    x = x(:);
    T = zeros(numel(x), deg + 1);
else
    T = casadi.SX.zeros(1, deg + 1);
    if isa(x, 'casadi.MX'), T = casadi.MX.zeros(1, deg + 1); end
end
T(:,1) = 1;
if deg >= 1, T(:,2) = x; end
for k = 3:deg + 1
    T(:,k) = 2*x.*T(:,k-1) - T(:,k-2);
end
end
