%% Reproduction of Malladi, Di Cairano, Weiss, "Nonlinear Model Predictive
%  Control of Coupled Rotational-Translational Spacecraft Relative Motion",
%  ACC 2019 (MERL TR2019-064).
%
%  Requires CasADi (with IPOPT) unzipped into external/casadi, e.g.
%  https://github.com/casadi/casadi/releases (casadi-*-osx-arm64-matlab2018b.zip).

projectRoot = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(projectRoot, 'src')));
addpath(fullfile(projectRoot, 'external', 'casadi'));
if ~exist('casadi.Opti', 'class')
    error('CasADi not found. Unzip the MATLAB release into external/casadi.');
end

p = paperParameters();
ctrl = buildNmpc(p);
log = runClosedLoop(p, ctrl);

resultsDir = fullfile(projectRoot, 'results');
if ~isfolder(resultsDir), mkdir(resultsDir); end
save(fullfile(resultsDir, 'closedLoop.mat'), 'log', 'p');
plotResults(log, p, resultsDir);
summarizeResults(log, p);
