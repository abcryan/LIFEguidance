%% Problem 4 phase sweep and per-spacecraft Problem 3 (run before main.m)
%  Sec. 5.2 "Additions for Problem 4": for fixed (tau = tf, phiTau) the
%  fleet problem separates into five single-spacecraft subproblems, eq.
%  (65). On a parallel pool, solve them on the 16-phase grid (and at
%  Problem 2's prescribed phase), each from two attitude-class guesses,
%  and solve Problem 3 for every spacecraft (J*_3,i, the price Pi4 of the
%  shared variables). main.m (Problem 4) starts the joint NLP from the
%  best grid points of this sweep. Writes results/problem4/sweep.mat and
%  sweep.json. About 10-15 min on 8 workers.
projectRoot = fileparts(mfilename('fullpath'));
addpath(projectRoot);
addpath(genpath(fullfile(projectRoot, 'src')));
setupCasadi();

p = lifeParameters();
ref = propagateReference(p, p.tf);
checkSunFeasibility(p, ref);            % same attitudes for all spacecraft
S = runFleetSweep(p, ref);

outDir = fullfile(projectRoot, 'results', 'problem4');
if ~isfolder(outDir), mkdir(outDir); end
save(fullfile(outDir, 'sweep.mat'), 'S', '-v7.3');
printFleetSweep(S, p);
writeSweepJson(S, p, fullfile(outDir, 'sweep.json'));

% Close the pool explicitly (MATLAB crashed shutting it down at exit)
delete(gcp('nocreate'));
