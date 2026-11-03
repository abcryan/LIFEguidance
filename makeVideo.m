%% Render the comparison animation of Problems 1-3 from a saved run of main.m
%  (and, optionally, the two per-problem animations of Sec. 5.3)
projectRoot = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(projectRoot, 'src')));

% Video settings
video.fps        = 30;                 % frames per second
video.duration   = 30;                 % [s] of video (constant playback speed)
video.height     = 1080;               % [px] frame height (1920x1080)
video.quality    = 95;                 % MPEG-4 quality, 0-100
video.perProblem = false;              % also the SER and formation videos per problem

resultsDir = fullfile(projectRoot, 'results');
pNow = lifeParameters();
R = {};
for k = 1:3                            % every saved problem, with its own p
    f = fullfile(resultsDir, sprintf('problem%d', k), 'solution.mat');
    if ~isfile(f), continue, end
    data = load(f, 'res', 'p');
    data.p.viz = pNow.viz;             % drawing settings always current
    data.res.p = data.p;
    R{end + 1} = data.res; %#ok<SAGROW>
end
p = pNow;                              % shared settings (Sun bounds, gates)
[~, iLong] = max(cellfun(@(r) r.sol.tau, R));
sched = videoSchedule(R{iLong}, 'fps', video.fps, 'duration', video.duration);
sched.height = video.height;
sched.quality = video.quality;
outDir = fullfile(resultsDir, 'comparison');
if ~isfolder(outDir), mkdir(outDir); end
videoCompare(R, p, sched, fullfile(outDir, 'video_comparison.mp4'));
fprintf('Video written to %s\n', outDir);

if video.perProblem
    for k = 1:numel(R)
        sk = videoSchedule(R{k}, 'fps', video.fps, 'duration', video.duration);
        sk.height = video.height; sk.quality = video.quality;
        dirK = fullfile(resultsDir, sprintf('problem%d', R{k}.problem));
        videoSER(R{k}, R{k}.p, sk, fullfile(dirK, 'video_SER.mp4'));
        videoFormation(R{k}, R{k}.p, sk, fullfile(dirK, 'video_formation.mp4'));
    end
end
