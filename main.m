%% Single-spacecraft guidance problems 1 and 2 at L2
%  Updated_LIFE_Guidance_Problem.pdf, Secs. 4-5: minimum-impulse acquisition
%  of a point on a rotating virtual formation near Sun-Earth L2 with 20
%  fixed throttleable thrusters. Problem 1: fixed phase, no path
%  constraint. Problem 2: adds the Sun-angle constraint (45). Direct Radau
%  collocation (CasADi/IPOPT) on one uniform mesh, independent forward
%  propagation, Sec. 5.3 figures per problem and a comparison.
%
%  CasADi (with IPOPT) is fetched into external/casadi by setupCasadi.
%  The video is rendered from the saved results by makeVideo.m.

problems = [1 2];                       % guidance problems to solve

projectRoot = fileparts(mfilename('fullpath'));
addpath(projectRoot);
addpath(genpath(fullfile(projectRoot, 'src')));
setupCasadi();

p = lifeParameters();
ref = propagateReference(p, p.tf);
resultsDir = fullfile(projectRoot, 'results');
for k = problems
    if k >= 2, checkSunFeasibility(p, ref); end
    ocp = buildOcp(p, ref, k);
    % Every problem starts from the same initial guess (Sec. 5.2), so the
    % solutions compare 1:1 (no warm start of Problem 2 from Problem 1)
    [sol, info] = solveOcp(ocp, initialGuess(ocp, p), p);
    res = validateSolution(sol, info, ref, p, k);
    outDir = fullfile(resultsDir, sprintf('problem%d', k));
    saveResults(res, ref, p, outDir);
    plotResults(res, p, outDir);
    printSummary(res, p);
end
compareProblems(resultsDir, p);
