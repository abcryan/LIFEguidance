function cmp = compareBurns(burns, prev, p)
%COMPAREBURNS Match burns of two refinement levels by time overlap and
%   compare their inertial force/torque integrals, duration and peak
%   thrust. cmp.ok is true when every burn of either level is matched and
%   all changes are within the p.mesh tolerances.
nb = numel(burns);
cmp.table = struct('burn', {}, 'match', {}, 'dF', {}, 'dT', {}, ...
                   'dDuration', {}, 'dPeak', {}, 'converged', {});
matchedPrev = false(1, numel(prev));
for b = 1:nb
    % Largest overlap; a near-impulsive burn may shift by about its own
    % duration between levels, so a gap up to the longer duration still matches
    best = 0; bestOv = -inf;
    for a = 1:numel(prev)
        ov = min(burns(b).t1, prev(a).t1) - max(burns(b).t0, prev(a).t0);
        if ov > bestOv && ov > -max(burns(b).duration, prev(a).duration)
            best = a; bestOv = ov;
        end
    end
    row.burn = b; row.match = best;
    if best == 0
        row.dF = inf; row.dT = inf; row.dDuration = inf; row.dPeak = inf;
    else
        A = prev(best); matchedPrev(best) = true;
        row.dF = norm(burns(b).Fc - A.Fc)/max(norm(A.Fc), p.mesh.burnFloorF);
        row.dT = norm(burns(b).Tc - A.Tc)/max(norm(A.Tc), p.mesh.burnFloorT);
        row.dDuration = abs(burns(b).duration - A.duration)/A.duration;
        row.dPeak = abs(burns(b).peakT - A.peakT)/A.peakT;
    end
    row.converged = row.dF < p.mesh.burnRelTol && row.dT < p.mesh.burnRelTol ...
        && row.dDuration < p.mesh.burnShapeTol && row.dPeak < p.mesh.burnShapeTol;
    cmp.table(b) = row;
end
cmp.unmatchedPrev = find(~matchedPrev);
cmp.ok = nb > 0 && all([cmp.table.converged]) && isempty(cmp.unmatchedPrev);
end
