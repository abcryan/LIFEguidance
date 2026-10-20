function R = alignZ(v)
%ALIGNZ Rotation taking e3 to the direction of v.
z = [0; 0; 1]; v = v/norm(v); c = z'*v;
if c < -1 + 1e-9, R = diag([1 -1 -1]); return, end
K = skew(cross(z, v));
R = eye(3) + K + K*K/(1 + c);
end
