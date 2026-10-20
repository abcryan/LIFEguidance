%% Weakly-constrained, fixed-phase, single-spacecraft guidance problem at L2
%  Guidance_Formulation_LIFE.pdf, Sec. 4: minimum-impulse acquisition of a
%  point on a rotating virtual formation near Sun-Earth L2 with 20 fixed
%  throttleable thrusters. Direct Radau collocation (CasADi/IPOPT), mesh
%  refinement, independent forward propagation and Sec. 4.4 outputs.
%
%  Stages can be switched off to reuse the saved results of a previous run.

projectRoot = fileparts(mfilename('fullpath'));
addpath(projectRoot);
addpath(genpath(fullfile(projectRoot, 'src')));
setupCasadi();

runSolve  = true;     % optimization + mesh refinement + validation (~2 min)
runPlots  = true;     % static figures
runVideos = true;     % two synchronized MPEG-4 animations (~10 min)

resultsDir = fullfile(projectRoot, 'results');
if ~isfolder(resultsDir), mkdir(resultsDir); end
resultsFile = fullfile(resultsDir, 'solution.mat');

p = lifeParameters();
if runSolve
    ref = propagateReference(p, p.tf);
    fprintf('Reference: Chebyshev fit error %.2e m on [0, tf]\n', ref.fitErr);
    [levels, final] = runMeshRefinement(p, ref, resultsDir);
    res = validateSolution(levels, final, ref, p);
    printSummary(res, p);
    saveResults(res, levels, ref, p, resultsDir);
else
    load(resultsFile, 'res', 'levels', 'p');
end

if runPlots
    plotResults(res, levels, p, resultsDir);
end
if runVideos
    sched = videoSchedule(res, p);
    videoSER(res, p, sched, fullfile(resultsDir, 'video_SER.mp4'));
    videoFormation(res, p, sched, fullfile(resultsDir, 'video_formation.mp4'));
end
