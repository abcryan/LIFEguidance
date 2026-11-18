function saveFleet(F, M, S, ref, p, outDir)
%SAVEFLEET Machine-readable outputs of Problem 4 (Sec. 5.3) in outDir:
%   solution.mat (fleet: per-spacecraft results with their own p, joint
%   solve, metrics; the sweep without its solution vectors; the reference
%   is in every scI/solution.mat), summary.json,
%   joint_starts.json, p3_per_spacecraft.json, and for every spacecraft
%   scI/ with the single-spacecraft outputs of saveResults and the figures
%   of plotResults (written in parallel, figures closed).
if ~isfolder(outDir), mkdir(outDir); end
nS = numel(F.res);
sweep = rmfield(S, 'sub');
fleet = F; metrics = M;
save(fullfile(outDir, 'solution.mat'), 'fleet', 'metrics', 'sweep', 'p', '-v7.3');

s.problem = 4;
s.status = F.info.status; s.iterations = F.info.iter; s.solve_time_s = F.info.wallTime;
s.joint_time_all_starts_s = F.jointTime; s.nlp_build_time_s = F.buildTime;
s.sweep_wall_time_s = S.wallTime; s.validation_time_s = F.validationTime;
s.N = p.ocp.N;
s.tau_s = M.tau; s.tau_at_upper_bound = F.info.tauAtUpper; s.tau_bound_multiplier = F.info.tauMultiplier;
s.phi_tau_deg = rad2deg(M.phiTau);
s.total_impulse_Ns = M.J4; s.propellant_kg = M.dm;
s.sum_J3i_Ns = sum(M.J3i); s.Pi4_Ns = M.Pi4; s.Pi4_relative = M.Pi4rel;
s.sweep_min_J4_Ns = M.sweepMin; s.sweep_min_phi_deg = rad2deg(M.sweepMinPhi);
s.prescribed_phase_bound_Ns = M.prescBound;
s.sensitivity_sum_Ns_per_rad = M.sensSum;
for i = nS:-1:1
    r = F.res{i};
    sc(i) = struct('spacecraft', i, 'name', p.fleet.name(i), 'J_Ns', M.Ji(i), ...
        'share', M.share(i), 'propellant_kg', M.dmi(i), 'J3_Ns', M.J3i(i), ...
        'phi3_deg', rad2deg(M.phi3i(i)), 'dJ_dphi_Ns_per_rad', M.dJdphi(i), ...
        'dJ_dphi_fd_Ns_per_rad', M.dJdphiFD(i), 'path_class', M.cls(i), ...
        'max_thrust_N', r.maxThrust, 'sun_max_violation_deg', r.sun.maxViolDeg, ...
        'sun_active_points', numel(r.sun.active), 'terminal', rmfield(r.term, {'pass', 'passAll'}), ...
        'terminal_pass', r.term.passAll, 'absolute_model_check_pass', r.absPass, ...
        'accepted', r.accepted);
end
s.spacecraft = sc;
s.min_distance = struct('m', M.dist.min, 'pair', M.dist.minPair, 't_s', M.dist.minTime, ...
    'pairs', M.dist.pairs, 'at_start_m', M.dist.atStart, 'at_arrival_m', M.dist.atArrival);
s.checks = M.checks;
s.expectations = M.expectations;
s.joint_starts = F.starts;
s.finite_differences = F.fd;
writeText(fullfile(outDir, 'summary.json'), jsonencode(s, 'PrettyPrint', true));
writeText(fullfile(outDir, 'joint_starts.json'), jsonencode(F.starts, 'PrettyPrint', true));
p3 = rmfield(S.p3, {'w', 'phiTau'});
for i = 1:nS, p3(i).phiTau_deg = rad2deg(S.p3(i).phiTau); end
writeText(fullfile(outDir, 'p3_per_spacecraft.json'), jsonencode(p3, 'PrettyPrint', true));

% Per spacecraft: the single-spacecraft outputs and figures
saveSpacecraft(F.res, ref, outDir);
end

function writeText(file, txt)
fid = fopen(file, 'w'); fprintf(fid, '%s', txt); fclose(fid);
end
