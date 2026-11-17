function key = physicsKey(p)
%PHYSICSKEY The spacecraft parameters compiled into the OCP functions
%   (ocpFunctions): mass, inertia, SRP coefficient, thrust bound and the
%   Sun-angle constraint. One compiled OCP serves every spacecraft with
%   the same key (setSpacecraft, buildFleetOcp).
key = [p.m; p.J(:); p.beta; p.Tmax; p.sun.n; p.sun.thetaMin; p.sun.thetaMax; p.sun.marginDeg];
end
