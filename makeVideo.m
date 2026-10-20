%% Render the two synchronized animations of Sec. 4.4 from a saved run of main.m
projectRoot = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(projectRoot, 'src')));

resultsDir = fullfile(projectRoot, 'results');
data = load(fullfile(resultsDir, 'solution.mat'), 'res', 'p');
sched = videoSchedule(data.res, data.p);
videoSER(data.res, data.p, sched, fullfile(resultsDir, 'video_SER.mp4'));
videoFormation(data.res, data.p, sched, fullfile(resultsDir, 'video_formation.mp4'));
fprintf('Videos written to %s\n', resultsDir);
