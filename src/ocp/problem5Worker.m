function out = problem5Worker(C, mode, varargin)
%PROBLEM5WORKER Problem 5 on a pool worker with the worker's compiled NLP
%   (parallel.pool.Constant C of buildFleetOcp with collision rows):
%   mode 'solve': problem5Solve(C.Value, varargin{:});
%   mode 'meta':  the row layout of the NLP (pairs, collision rows) and
%                 its size, which also builds it on this worker.
fo = C.Value;
switch mode
    case 'solve'
        out = problem5Solve(fo, varargin{:});
    case 'meta'
        out = struct('pairs', fo.pairs, 'iCol', fo.iCol, 'colPair', fo.colPair, ...
                     'colTime', fo.colTime, 'nW', fo.nW, 'nG', fo.nG);
end
end
