function [F, M] = runProblem5(p, ref, resultsDir)
%RUNPROBLEM5 Problem 5 in main.m: loads the Problem 4 solution
%   (results/problem4/solution.mat, which must match the current
%   settings), checks the formations (step (i)), solves Problem 5
%   (solveProblem5), verifies it (problem5Metrics) and writes the outputs
%   and figures to results/problem5/.
f = fullfile(resultsDir, 'problem4', 'solution.mat');
if ~isfile(f)
    error('LIFE:noProblem4', 'Problem 5 starts from Problem 4: run main with problem 4 first (%s missing).', f);
end
P4 = load(f, 'fleet', 'metrics', 'p');
if ~strcmp(fleetKey(P4.p), fleetKey(p))
    error('LIFE:staleProblem4', '%s was computed with other settings than lifeParameters(); rerun Problem 4.', f);
end
tAll = tic;
feas = checkCollisionFeasibility(p);
F = solveProblem5(p, ref, P4);
M = problem5Metrics(F, P4, p, feas);
outDir = fullfile(resultsDir, 'problem5');
saveProblem5(F, M, ref, p, outDir);
plotFleet(F, struct('tau', M.tau, 'cls', M.cls), [], p, outDir, ...
          {'fleet_trajectories', 'fleet_tracking', 'fleet_sun', 'fleet_thrust'});
plotProblem5(F, M, P4, p, outDir);
printProblem5(F, M, p);
fprintf('Problem 5 (starts, sweep, checks, validation, outputs): %.1f min\n', toc(tAll)/60);
end
