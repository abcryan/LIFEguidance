function q = dcmToQuat(R)
%DCMTOQUAT Unit quaternion (scalar first, q(1) >= 0) of a rotation matrix.
%   Shepperd's method.
tr = trace(R);
[~, i] = max([tr, R(1,1), R(2,2), R(3,3)]);
switch i
    case 1
        w = 0.5*sqrt(1 + tr);
        q = [w; (R(3,2) - R(2,3))/(4*w); (R(1,3) - R(3,1))/(4*w); (R(2,1) - R(1,2))/(4*w)];
    case 2
        x = 0.5*sqrt(1 + 2*R(1,1) - tr);
        q = [(R(3,2) - R(2,3))/(4*x); x; (R(1,2) + R(2,1))/(4*x); (R(1,3) + R(3,1))/(4*x)];
    case 3
        y = 0.5*sqrt(1 + 2*R(2,2) - tr);
        q = [(R(1,3) - R(3,1))/(4*y); (R(1,2) + R(2,1))/(4*y); y; (R(2,3) + R(3,2))/(4*y)];
    case 4
        z = 0.5*sqrt(1 + 2*R(3,3) - tr);
        q = [(R(2,1) - R(1,2))/(4*z); (R(1,3) + R(3,1))/(4*z); (R(2,3) + R(3,2))/(4*z); z];
end
if q(1) < 0, q = -q; end
q = q/norm(q);
end
