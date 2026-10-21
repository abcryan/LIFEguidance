%% Render the two synchronized animations of Sec. 4.4 from a saved run of main.m
projectRoot = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(projectRoot, 'src')));

% Video settings
video.fps      = 30;                   % frames per second
video.duration = 30;                   % [s] of video (constant playback speed)
video.height   = 1080;                 % [px] frame height (1920x1080)
video.quality  = 95;                   % MPEG-4 quality, 0-100

resultsDir = fullfile(projectRoot, 'results');
data = load(fullfile(resultsDir, 'solution.mat'), 'res', 'p');
p = data.p;
pNow = lifeParameters();
p.viz = pNow.viz;                      % drawing settings always current
sched = videoSchedule(data.res, 'fps', video.fps, 'duration', video.duration);
sched.height = video.height;
sched.quality = video.quality;
videoSER(data.res, p, sched, fullfile(resultsDir, 'video_SER.mp4'));
videoFormation(data.res, p, sched, fullfile(resultsDir, 'video_formation.mp4'));
fprintf('Videos written to %s\n', resultsDir);
