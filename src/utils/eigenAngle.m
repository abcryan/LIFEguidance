function [thAcos, thAtan] = eigenAngle(q)
%EIGENANGLE Eigenaxis angle of the rotation represented by quaternion q.
%   thAcos is the form of eqs. (60) and (5), evaluated on R(q) without
%   renormalizing q (argument clamped to [-1, 1]); thAtan =
%   2*atan2(|q_v|, |q_w|) has no ~1.5e-8 rad round-off floor.
%   q may be 4xn.
n = size(q, 2);
thAcos = zeros(1, n); thAtan = zeros(1, n);
for j = 1:n
    w = q(1,j); v = q(2:4,j);
    trR = 3*w^2 - v'*v;                 % trace of R(q), no normalization
    thAcos(j) = acos(min(max((trR - 1)/2, -1), 1));
    thAtan(j) = 2*atan2(norm(v), abs(w));
end
end
