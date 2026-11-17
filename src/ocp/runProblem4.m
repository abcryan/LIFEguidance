function [F, M] = runProblem4(p, ref, resultsDir)
%RUNPROBLEM4 Problem 4 in main.m: loads the sweep of fleetSweep.m
%   (results/problem4/sweep.mat, which must match the current settings),
%   solves the joint NLP, checks and validates it (solveFleet), and
%   writes the fleet outputs and figures to results/problem4/.
outDir = fullfile(resultsDir, 'problem4');
f = fullfile(outDir, 'sweep.mat');
if ~isfile(f)
    error('LIFE:noSweep', 'Problem 4 needs the sweep: run fleetSweep.m first (%s missing).', f);
end
d = load(f, 'S');
if ~strcmp(d.S.key, fleetKey(p))
    error('LIFE:staleSweep', ['%s was computed with other settings than lifeParameters(); ' ...
          'rerun fleetSweep.m.'], f);
end
tAll = tic;
F = solveFleet(p, ref, d.S);
M = fleetMetrics(F, d.S, p);
saveFleet(F, M, d.S, ref, p, outDir);
plotFleet(F, M, d.S, p, outDir);
printFleetSummary(F, M, p);
fprintf('Problem 4 (joint solve, checks, validation, outputs): %.1f min\n', toc(tAll)/60);
end
