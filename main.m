%% Weakly-constrained, fixed-phase, single-spacecraft guidance problem at L2
%  Guidance_Formulation_LIFE.pdf, Sec. 4: minimum-impulse acquisition of a
%  point on a rotating virtual formation near Sun-Earth L2 with 20 fixed
%  throttleable thrusters. Direct Radau collocation (CasADi/IPOPT) on one
%  uniform mesh, independent forward propagation, Sec. 4.4 figures.
%
%  CasADi (with IPOPT) is fetched into external/casadi by setupCasadi.
%  The videos are rendered from the saved results by makeVideo.m.

projectRoot = fileparts(mfilename('fullpath'));
addpath(projectRoot);
addpath(genpath(fullfile(projectRoot, 'src')));
setupCasadi();

p = lifeParameters();
ref = propagateReference(p, p.tf);
ocp = buildOcp(p, ref);
[sol, info] = solveOcp(ocp, initialGuess(ocp, p), p);
res = validateSolution(sol, info, ref, p);

resultsDir = fullfile(projectRoot, 'results');
if ~isfolder(resultsDir), mkdir(resultsDir); end
saveResults(res, ref, p, resultsDir);
plotResults(res, p, resultsDir);
printSummary(res, p);
