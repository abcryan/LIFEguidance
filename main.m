%% Single-spacecraft guidance problems 1, 2 and 3 at L2
%  Updated LIFE guidance problem, Secs. 4-5: minimum-impulse acquisition of
%  a point on a rotating virtual formation near Sun-Earth L2 with 20 fixed
%  throttleable thrusters. Problem 1: fixed phase, no path constraint.
%  Problem 2: adds the Sun-angle constraint (45). Problem 3: adds the free
%  arrival phase phiTau, solved by multi-start. Direct Radau collocation
%  (CasADi/IPOPT) on one uniform mesh, independent forward propagation,
%  Sec. 5.3 figures per problem and a comparison.
%
%  CasADi (with IPOPT) is fetched into external/casadi by setupCasadi.
%  The video is rendered from the saved results by makeVideo.m, the
%  Problem 3 phase sweep by phaseSweep.m.

problems = [1 2 3];                     % guidance problems to solve

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
    % solutions compare 1:1; Problem 3 from every start phase
    starts = [];
    if k < 3
        [sol, info] = solveOcp(ocp, initialGuess(ocp, p), p);
    else
        [sol, info, starts] = solveMultiStart(ocp, p);
    end
    pk = p;                             % the solution's formation phase
    pk.target.phase0 = sol.phase0; pk.phi0 = sol.phase0;
    res = validateSolution(sol, info, ref, pk, k);
    res.multistart = starts;
    outDir = fullfile(resultsDir, sprintf('problem%d', k));
    saveResults(res, ref, pk, outDir);
    plotResults(res, pk, outDir);
    printSummary(res, pk);
end
compareProblems(resultsDir, p);
