%% Render the comparison animation of Problems 1-3 and the Problem 4 fleet
%  animation from a saved run of main.m (and, optionally, the two
%  per-problem animations of Sec. 5.3)
projectRoot = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(projectRoot, 'src')));

% Video settings
video.fps        = 29;                 % frames per second
video.duration   = 30;                 % [s] of video (constant playback speed)
video.height     = 1080;               % [px] frame height (1920x1080)
video.quality    = 100;                % MPEG-4 quality, 0-100
video.perProblem = false;              % also the SER and formation videos per problem
video.compare    = true;               % Problems 1-3: video_comparison.mp4
video.fleet      = true;               % Problem 4: problem4/video_fleet.mp4

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
if video.compare && ~isempty(R)
    [~, iLong] = max(cellfun(@(r) r.sol.tau, R));
    sched = videoSchedule(R{iLong}, 'fps', video.fps, 'duration', video.duration);
    sched.height = video.height;
    sched.quality = video.quality;
    outDir = fullfile(resultsDir, 'comparison');
    if ~isfolder(outDir), mkdir(outDir); end
    videoCompare(R, p, sched, fullfile(outDir, 'video_comparison.mp4'));
    fprintf('Video written to %s\n', outDir);
end

f4 = fullfile(resultsDir, 'problem4', 'solution.mat');
if video.fleet && isfile(f4)
    d = load(f4, 'fleet');
    for i = 1:numel(d.fleet.res), d.fleet.res{i}.p.viz = pNow.viz; end
    sk = videoSchedule(d.fleet.res{1}, 'fps', video.fps, 'duration', video.duration);
    sk.height = video.height; sk.quality = video.quality;
    videoFleet(d.fleet, p, sk, fullfile(resultsDir, 'problem4', 'video_fleet.mp4'));
    fprintf('Video written to %s\n', fullfile(resultsDir, 'problem4'));
end

if video.perProblem
    for k = 1:numel(R)
        sk = videoSchedule(R{k}, 'fps', video.fps, 'duration', video.duration);
        sk.height = video.height; sk.quality = video.quality;
        dirK = fullfile(resultsDir, sprintf('problem%d', R{k}.problem));
        videoSER(R{k}, R{k}.p, sk, fullfile(dirK, 'video_SER.mp4'));
        videoFormation(R{k}, R{k}.p, sk, fullfile(dirK, 'video_formation.mp4'));
    end
end
