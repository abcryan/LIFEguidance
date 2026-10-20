function r = chebEvalRef(ref, t)
%CHEBEVALREF Fitted reference position r^b_bf(t) [m]. Double t (1xn) or
%   CasADi scalar t.
x = 2*t/ref.tEnd - 1;
if isa(t, 'double')
    r = ref.r0 + (chebBasis(x(:), ref.deg)*ref.coef)';
else
    r = ref.r0 + (chebBasis(x, ref.deg)*ref.coef)';
end
end
