function theta = eigenaxisAngle(R)
%EIGENAXISANGLE Relative eigenaxis attitude error, eq. (5) [rad].
theta = acos(min(max((trace(R) - 1)/2, -1), 1));
end
