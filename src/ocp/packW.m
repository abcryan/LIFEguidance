function w = packW(ocp, tauH, X, XC, U)
%PACKW Decision vector from SI states X (13xN+1), XC (13xNd), T/Tmax U.
w = zeros(ocp.nW, 1);
w(ocp.ix.tau) = tauH;
w(ocp.ix.X)  = reshape(X./ocp.S, [], 1);
w(ocp.ix.XC) = reshape(XC./ocp.S, [], 1);
w(ocp.ix.U)  = U(:);
end
