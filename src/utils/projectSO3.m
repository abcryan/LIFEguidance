function R = projectSO3(R)
%PROJECTSO3 Closest rotation matrix (removes integration drift).
[U,~,V] = svd(R);
R = U*diag([1 1 det(U*V')])*V';
end
