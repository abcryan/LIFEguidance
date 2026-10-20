function [f, tau] = thrusterWrench(T, p)
%THRUSTERWRENCH Body-frame force and torque, eqs. (29)-(30). T is 20xn.
f = p.Bf*T;
tau = p.Btau*T;
end
