function w = packW(ocp, tauH, X, XC, U, phiTau, w)
%PACKW Decision vector from SI states X (13xN+1), XC (13xNd), T/Tmax U and,
%   for Problems 3-4, the arrival phase phiTau [rad]. With w given, the
%   block of ocp (e.g. one spacecraft view of the fleet OCP) is written
%   into it; otherwise into zeros.
if nargin < 7 || isempty(w), w = zeros(ocp.nW, 1); end
w(ocp.ix.tau) = tauH;
w(ocp.ix.X)  = reshape(X./ocp.S, [], 1);
w(ocp.ix.XC) = reshape(XC./ocp.S, [], 1);
w(ocp.ix.U)  = U(:);
if ~isempty(ocp.ix.phi), w(ocp.ix.phi) = phiTau/pi; end
end
