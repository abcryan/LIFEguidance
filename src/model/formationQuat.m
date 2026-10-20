function q = formationQuat(t, p)
%FORMATIONQUAT Quaternion of the target formation attitude R^c_f(t),
%   continuous in t.
[~, q] = formationAttitude(p.target, t);
end
