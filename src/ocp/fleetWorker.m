function r = fleetWorker(C, varargin)
%FLEETWORKER solveSubproblem on a pool worker with the worker's compiled
%   parametric Problem 3 OCP (parallel.pool.Constant C).
r = solveSubproblem(C.Value, varargin{:});
end
