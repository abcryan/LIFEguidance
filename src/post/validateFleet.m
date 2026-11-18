function [res, tV] = validateFleet(sols, info, ref, p, problem)
%VALIDATEFLEET Validation of every spacecraft of a joint fleet solution
%   (Problems 4 and 5) with the shared formation phase: validateSolution
%   (forward propagation, arrival gates, Sun angle, absolute model, all
%   unchanged) in parallel on the pool, res{i} with its own p and index.
%   sols{i}: single-spacecraft solutions (solveFleetOcp); info: the joint
%   solve's info (with x, lamX, lamG). tV: wall time [s].
nS = numel(sols);
pk = cell(1, nS); inf_i = cell(1, nS);
for i = 1:nS
    pk{i} = spacecraftParams(p, i);
    pk{i}.target.phase0 = sols{i}.phase0; pk{i}.phi0 = sols{i}.phase0;
    in = rmfield(info, intersect(fieldnames(info), {'x', 'lamX', 'lamG'}));
    in.boundViolation = max([0; -sols{i}.U(:); sols{i}.U(:) - 1])*p.Tmax;
    inf_i{i} = in;
end
fleetPool(p);
tV = tic;
res = cell(1, nS);
parfor i = 1:nS
    res{i} = validateSolution(sols{i}, inf_i{i}, ref, pk{i}, problem);
end
for i = 1:nS
    res{i}.p = pk{i};
    res{i}.sc = i;
end
tV = toc(tV);
end
