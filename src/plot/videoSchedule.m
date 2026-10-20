function sched = videoSchedule(res, opts)
%VIDEOSCHEDULE Uniform frame times shared by both videos: constant
%   playback speed-up tau/duration over [0, tau].
arguments
    res struct
    opts.fps (1,1) double = 15
    opts.duration (1,1) double = 30      % [s] of video
end
n = round(opts.duration*opts.fps) + 1;
sched.t = linspace(0, res.sol.tau, n);
sched.fps = opts.fps;
sched.speed = res.sol.tau/opts.duration;
sched.label = sprintf('playback x%.0f', sched.speed);
end
