function a = gravity(r, mu)
%GRAVITY Two-body point-mass gravitational acceleration.
a = -mu*r/(sqrt(r'*r))^3;
end
