function T = thrustInWindow(res, t0, t1)
%THRUSTINWINDOW Peak thrust of each thruster over [t0, t1] (20x1), so that
%   a burn shorter than one video frame step still shows up in its frame.
tn = res.sol.tNode;
k = find(tn(1:end-1) < t1 & tn(2:end) > t0);
if isempty(k), k = min(max(discretize(t1, tn), 1), res.sol.N); end
T = max(res.sol.T(:, k), [], 2);
end
