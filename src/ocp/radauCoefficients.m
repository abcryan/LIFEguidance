function [cRoot, C, D] = radauCoefficients(d)
%RADAUCOEFFICIENTS Lagrange collocation data on [0, 1] with d Radau points.
%   cRoot = [0, c_1..c_d]; C(j,r) = dL_j/ds at cRoot(r); D(j) = L_j(1).
cRoot = [0, casadi.collocation_points(d, 'radau')];
C = zeros(d + 1); D = zeros(d + 1, 1);
for j = 1:d + 1
    L = 1;
    for r = [1:j-1, j+1:d+1]
        L = conv(L, [1, -cRoot(r)])/(cRoot(j) - cRoot(r));
    end
    D(j) = polyval(L, 1);
    dL = polyder(L);
    for r = 1:d + 1
        C(j,r) = polyval(dL, cRoot(r));
    end
end
end
