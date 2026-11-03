function w = packW(ocp, tauH, X, XC, U, phiTau)
%PACKW Decision vector from SI states X (13xN+1), XC (13xNd), T/Tmax U and,
%   for Problem 3, the arrival phase phiTau [rad].
w = zeros(ocp.nW, 1);
w(ocp.ix.tau) = tauH;
w(ocp.ix.X)  = reshape(X./ocp.S, [], 1);
w(ocp.ix.XC) = reshape(XC./ocp.S, [], 1);
w(ocp.ix.U)  = U(:);
if ~isempty(ocp.ix.phi), w(ocp.ix.phi) = phiTau/pi; end
end
