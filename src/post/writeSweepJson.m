function writeSweepJson(S, p, file)
%WRITESWEEPJSON The Problem 4 sweep without solution vectors, as JSON.
nS = size(S.Jbest, 2);
kept = ["guess via +z_c", "guess via -z_c", "retry, direct slerp", "retry, warm from nearest phase"];
for m = numel(S.phis):-1:1
    for i = nS:-1:1
        sub(i) = struct('spacecraft', i, 'J_Ns', S.Jbest(m, i), 'class', S.clsBest(m, i), ...
            'kept_solve', kept(S.vBest(m, i)), 'dJ_dphi_Ns_per_rad', S.dJbest(m, i), ...
            'J_both_guesses_Ns', squeeze(S.J(m, i, :))', ...
            'status_both_guesses', squeeze(S.status(m, i, :))', ...
            'class_both_guesses', squeeze(S.cls(m, i, :))', ...
            'iterations_both_guesses', squeeze(S.iter(m, i, :))');
    end
    rows(m) = struct('phiTau_deg', rad2deg(S.phis(m)), 'prescribed_P2', S.isPresc(m), ...
        'J4_Ns', S.J4(m), 'joint_start', any(S.gridBest == m), 'spacecraft', sub);
end
p3 = rmfield(S.p3, {'w', 'phiTau'});
for i = 1:nS
    p3(i).phiTau_deg = rad2deg(S.p3(i).phiTau);
    p3(i).name = p.fleet.name(i);
end
out = struct('sweep', rows, 'retries', S.retry, 'problem3_per_spacecraft', p3, ...
    'wall_time_s', S.wallTime, 'solver_time_s', S.cpuTime);
fid = fopen(file, 'w'); fprintf(fid, '%s', jsonencode(out, 'PrettyPrint', true)); fclose(fid);
end
