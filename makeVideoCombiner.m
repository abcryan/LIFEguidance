%% Render the per-problem animations of the combiner's Problem 3
%  From the saved run of problem3Combiner.m: video_SER.mp4 and
%  video_formation.mp4 in results/problem3_combiner/, with the settings
%  of makeVideo.m.
projectRoot = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(projectRoot, 'src')));

video.fps      = 29;                   % frames per second
video.duration = 30;                   % [s] of video (constant playback speed)
video.height   = 1080;                 % [px] frame height (1920x1080)
video.quality  = 100;                  % MPEG-4 quality, 0-100

outDir = fullfile(projectRoot, 'results', 'problem3_combiner');
f = fullfile(outDir, 'solution.mat');
if ~isfile(f)
    error('LIFE:noCombiner', 'Run problem3Combiner.m first (%s missing).', f);
end
d = load(f, 'res', 'p');
d.p.viz = lifeParameters().viz;        % drawing settings always current
d.res.p = d.p;
sk = videoSchedule(d.res, 'fps', video.fps, 'duration', video.duration);
sk.height = video.height; sk.quality = video.quality;
videoSER(d.res, d.p, sk, fullfile(outDir, 'video_SER.mp4'));
videoFormation(d.res, d.p, sk, fullfile(outDir, 'video_formation.mp4'));
fprintf('Videos written to %s\n', outDir);
