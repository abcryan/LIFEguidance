%% Reproduction of Malladi, Di Cairano, Weiss, "Nonlinear Model Predictive
%  Control of Coupled Rotational-Translational Spacecraft Relative Motion",
%  ACC 2019 (MERL TR2019-064).
%
%  CasADi (with IPOPT) is fetched into external/casadi by setupCasadi.

projectRoot = fileparts(mfilename('fullpath'));
addpath(projectRoot);
addpath(genpath(fullfile(projectRoot, 'src')));
setupCasadi();

p = paperParameters();
ctrl = buildNmpc(p);
log = runClosedLoop(p, ctrl);

resultsDir = fullfile(projectRoot, 'results');
if ~isfolder(resultsDir), mkdir(resultsDir); end
save(fullfile(resultsDir, 'closedLoop.mat'), 'log', 'p');
plotResults(log, p, resultsDir);
summarizeResults(log, p);
