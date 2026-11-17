function pool = fleetPool(p)
%FLEETPOOL Process-based parallel pool with p.par.nWorkers workers (the
%   current pool if one is open), with src/ and CasADi on every worker's
%   path. Used by the Problem 4 sweep and validation.
pool = gcp('nocreate');
if isempty(pool)
    pool = parpool('Processes', p.par.nWorkers);
end
% No idle shutdown: the client may compute for a long time between two
% parallel stages (the default 30 min timeout closed the pool once)
pool.IdleTimeout = Inf;
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
wait(parfevalOnAll(pool, @addpath, 0, genpath(fullfile(root, 'src')), ...
                   fullfile(root, 'external', 'casadi')));
end
