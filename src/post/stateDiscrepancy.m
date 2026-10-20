function d = stateDiscrepancy(xa, xb)
%STATEDISCREPANCY Pointwise differences of two relative state histories
%   (13xn): position, velocity, attitude (angle of R_a' R_b) and rate.
d.r = vecnorm(xa(8:10,:) - xb(8:10,:));
d.v = vecnorm(xa(11:13,:) - xb(11:13,:));
d.w = vecnorm(xa(5:7,:) - xb(5:7,:));
n = size(xa, 2);
d.theta = zeros(1, n);
for k = 1:n
    qa = xa(1:4,k)/norm(xa(1:4,k)); qb = xb(1:4,k)/norm(xb(1:4,k));
    dq = quatMul([qa(1); -qa(2:4)], qb);
    [~, d.theta(k)] = eigenAngle(dq);
end
d.max = struct('r', max(d.r), 'v', max(d.v), 'theta', max(d.theta), 'w', max(d.w));
end
