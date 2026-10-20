function R = quatToDcm(q)
%QUATTODCM Rotation matrix of a Hamilton, scalar-first quaternion.
%   R = R(q)/|q|^2 is exactly orthogonal for any nonzero q, so a drift of
%   |q| never distorts R; the drift itself is monitored separately.
%   Double or CasADi.
w = q(1); v = q(2:4);
R = ((w^2 - v'*v)*eye(3) + 2*(v*v') + 2*w*skew(v))/(q'*q);
end
