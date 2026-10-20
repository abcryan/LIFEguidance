function dz = relativeRhs(t, z, T, p)
%RELATIVERHS Scaled relative state (13) plus normalized reference (6), as
%   integrated by the independent propagation.
S = p.ocp.xScale;
xd = relativeDynamics(t, S.*z(1:13), T, z(14:16)*p.dse, p);
dz = [xd./S; p.Om*cr3bpNormalized(z(14:19), p.eta)];
end
