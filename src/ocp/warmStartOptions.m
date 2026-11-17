function o = warmStartOptions()
%WARMSTARTOPTIONS IPOPT options for a primal-dual warm start from a nearby
%   optimum (Problem 4: the per-spacecraft Problem 3 from its best sweep
%   point, fleetSweep.m; for moving the phase, see continuationOptions).
%   With IPOPT's default initialization, the bound push moves the start
%   off the optimum and the solve may end in another local optimum.
o = struct('warm_start_init_point', 'yes', 'warm_start_bound_push', 1e-9, ...
    'warm_start_bound_frac', 1e-9, 'warm_start_slack_bound_push', 1e-9, ...
    'warm_start_slack_bound_frac', 1e-9, 'warm_start_mult_bound_push', 1e-9, 'mu_init', 1e-6);
end
