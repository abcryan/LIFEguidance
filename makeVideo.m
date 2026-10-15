%% Render the inertial-frame animation (cf. Fig. 5) from a saved run of main.m
projectRoot = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(projectRoot, 'src')));

data = load(fullfile(projectRoot, 'results', 'closedLoop.mat'));
videoFile = fullfile(projectRoot, 'results', 'inertialFrame.mp4');
animateInertial(data.log, data.p, videoFile);
fprintf('Video written to %s\n', videoFile);
