function key = fleetKey(p)
%FLEETKEY Settings that determine the Problem 4 sweep (fleetSweep.m), as a
%   string; main.m reuses a saved sweep only if its key is equal.
k.N = p.ocp.N; k.d = p.ocp.d; k.Tmax = p.Tmax; k.tf = p.tf; k.tauMin = p.tauMin;
k.target = rmfield(p.target, {'Rstar', 'qstar'});
k.init = rmfield(p.init, {'Rstar', 'qstar'});
k.fleet = p.fleet; k.sun = p.sun; k.m = p.m; k.J = p.J; k.beta = p.beta;
k.haloIC = p.haloIC; k.xScale = p.ocp.xScale; k.Jscale = p.ocp.Jscale;
k.ipopt = rmfield(p.ocp.ipopt, {'print_level', 'print_frequency_iter'});
k.p4 = rmfield(p.ocp.p4, intersect(fieldnames(p.ocp.p4), ...
    {'fdStepDeg'}));                     % not used by the sweep
k.p3Starts = p.ocp.p3Starts; k.uGuess = p.ocp.uGuess;
key = jsonencode(k);
end
