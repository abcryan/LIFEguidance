%% Render the two synchronized animations of Sec. 4.4 from a saved run of main.m
projectRoot = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(projectRoot, 'src')));

resultsDir = fullfile(projectRoot, 'results');
data = load(fullfile(resultsDir, 'solution.mat'), 'res', 'p');
p = data.p;
pNow = lifeParameters();
p.viz = pNow.viz;                      % drawing settings always current
sched = videoSchedule(data.res);       % uniform playback speed
videoSER(data.res, p, sched, fullfile(resultsDir, 'video_SER.mp4'));
videoFormation(data.res, p, sched, fullfile(resultsDir, 'video_formation.mp4'));
fprintf('Videos written to %s\n', resultsDir);
