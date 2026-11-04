function sigma = sunCosine(t, q, rbi, p)
%SUNCOSINE Cosine of the Sun angle, eq. (47): sigma = rhat_is' R^b_i n^i,
%   with R^b_i = R^b_f(t) R(q) and rhat_is the unit vector from the
%   position rbi (SER, [m]) to the Sun. The NLP passes the formation centre
%   r_bf for rbi (Sec. 5.2, error ~1e-9 rad), the verification the
%   spacecraft position. q may be 4xn with t 1xn and rbi 3xn (double), or
%   a single CasADi sample.
n = size(q, 2);
sigma = q(1,:);                         % right type and size
for k = 1:n
    fr = frames(t(k), p);
    s = p.rbs - rbi(:,k);
    sigma(k) = (s'/sqrt(s'*s))*(fr.Rbf*(quatToDcm(q(:,k))*p.sun.n));
end
end
