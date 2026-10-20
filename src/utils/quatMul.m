function r = quatMul(p, q)
%QUATMUL Hamilton product p (x) q, scalar first. Double or CasADi.
r = [p(1)*q(1) - p(2:4)'*q(2:4);
     p(1)*q(2:4) + q(1)*p(2:4) + cross(p(2:4), q(2:4))];
end
