function sched = videoSchedule(res, p, opts)
%VIDEOSCHEDULE Event-aware frame times shared by both videos: coasts at a
%   high constant speed-up, every burn (padded) in slow motion so that it
%   spans >= opts.burnVideo seconds of video, geometric speed ramps of
%   ~opts.rampVideo seconds and a short hold at each burn start.
arguments
    res struct
    p struct
    opts.fps (1,1) double = 15
    opts.coastVideo (1,1) double = 12    % [s] of video for the whole span
    opts.burnVideo (1,1) double = 1
    opts.rampVideo (1,1) double = 0.3
    opts.holdVideo (1,1) double = 0.2
end
tau = res.sol.tau;
vCoast = tau/opts.coastVideo;
burns = res.burns;
W = zeros(0, 3);                         % [t0 t1 speed] slow-motion windows
for b = 1:numel(burns)
    d = burns(b).duration;
    v = min(vCoast, max(d/opts.burnVideo, 0.25));
    W(end+1,:) = [max(burns(b).t0 - 0.5*d, 0), min(burns(b).t1 + 0.5*d, tau), v]; %#ok<AGROW>
end
speed = @(t) speedAt(t, W, vCoast, opts.rampVideo);

t = 0; ts = [];
while t < tau
    ts(end+1) = t; %#ok<AGROW>
    t = t + speed(t)/opts.fps;
end
ts(end+1) = tau;
sp = arrayfun(speed, ts);

% Holds at burn starts
tAll = []; spAll = []; hold = [];
nHold = round(opts.holdVideo*opts.fps);
for k = 1:numel(ts)
    tAll(end+1) = ts(k); spAll(end+1) = sp(k); hold(end+1) = 0; %#ok<AGROW>
    for b = 1:numel(burns)
        if ts(k) <= burns(b).t0 && (k == numel(ts) || ts(k+1) > burns(b).t0)
            tAll(end+(1:nHold)) = burns(b).t0;
            spAll(end+(1:nHold)) = 0;
            hold(end+(1:nHold)) = b;
        end
    end
end
sched.t = tAll; sched.speed = spAll; sched.hold = hold;
sched.fps = opts.fps; sched.vCoast = vCoast; sched.windows = W;
sched.label = strings(1, numel(tAll));
for k = 1:numel(tAll)
    inBurn = find(tAll(k) >= [burns.t0] & tAll(k) <= [burns.t1], 1);
    if hold(k) > 0
        sched.label(k) = sprintf('paused: burn %d starts', hold(k));
    elseif spAll(k) < 0.99*vCoast
        if isempty(inBurn)
            sched.label(k) = sprintf('x%.3g slow motion', spAll(k));
        else
            sched.label(k) = sprintf('x%.3g slow motion: burn %d', spAll(k), inBurn);
        end
    else
        sched.label(k) = sprintf('x%.0f', spAll(k));
    end
end
end

function v = speedAt(t, W, vCoast, rampVideo)
v = vCoast;
for i = 1:size(W, 1)
    vw = W(i,3);
    rampPhys = rampVideo*sqrt(vw*vCoast);
    dist = max([W(i,1) - t, t - W(i,2), 0]);
    if dist == 0
        vi = vw;
    elseif dist < rampPhys
        vi = vw*(vCoast/vw)^(dist/rampPhys);
    else
        vi = vCoast;
    end
    v = min(v, vi);
end
end
