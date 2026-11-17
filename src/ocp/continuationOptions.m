function o = continuationOptions()
%CONTINUATIONOPTIONS IPOPT options for a primal start from a nearby
%   solution that must stay in its local optimum while the shared phase
%   moves by a few degrees (Problem 4: the joint NLP from the sweep grid
%   points, the finite-difference check): the start is kept (bound pushes
%   1e-8 instead of 1e-2, which lifts every thruster off zero) and the
%   barrier starts small (mu_init 1e-5), with IPOPT's own dual
%   initialization. Measured at N = 400, d = 5 on the five subproblems
%   stepped by 0.1 and 2 deg from the 0 deg sweep solution: all 10 solved
%   in 6-43 s and stayed in their local optimum (cost within the first-order
%   prediction). The alternatives fail: default initialization changes
%   the local optimum after 0.1 deg (collector 4 145.37 -> 147.02 N s) and
%   ends the joint NLP in worse optima (1405.9 and 1423.2 N s from grid
%   points at 1393.3 and 1386.1 N s); a primal-dual warm start
%   (warmStartOptions) needs 140-260 iterations for 0.1 deg and stalls for
%   2 deg, and the joint NLP from it did not converge in 900 s.
o = struct('bound_push', 1e-8, 'bound_frac', 1e-8, 'slack_bound_push', 1e-8, ...
    'slack_bound_frac', 1e-8, 'mu_init', 1e-5);
end
