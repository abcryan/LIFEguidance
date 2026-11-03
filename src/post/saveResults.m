function saveResults(res, ref, p, outDir)
%SAVERESULTS Machine-readable outputs for reproducibility (Sec. 5.3).
%   solution.mat (v7.3/HDF5), summary.json, settings.json and CSV files.
if ~isfolder(outDir), mkdir(outDir); end
ref = rmfield(ref, 'eval');             % function handles are not portable
save(fullfile(outDir, 'solution.mat'), 'res', 'ref', 'p', '-v7.3');

s.problem = res.problem;
s.status = res.info.status;
s.iterations = res.info.iter;
s.solve_time_s = res.info.wallTime;
s.accepted = res.accepted;
s.N = res.sol.N;
s.tau_s = res.sol.tau;
s.tau_at_upper_bound = res.info.tauAtUpper;
s.tau_bound_multiplier = res.info.tauMultiplier;
s.phi_tau_deg = rad2deg(res.sol.phiTau);
s.phi_tau_free = res.problem >= 3;
s.phase0_deg = rad2deg(res.sol.phase0);
s.total_impulse_Ns = res.J;
s.propellant_kg = res.dm;
s.max_thrust_N = res.maxThrust;
s.thrust_bound_violation_N = res.boundViolation;
s.terminal = rmfield(res.term, {'pass', 'passAll'});
s.terminal_pass = res.term.pass;
s.terminal_pass_all = res.term.passAll;
s.prop_vs_collocation = res.disc.max;
s.absolute_model_check = res.absDisc.max;
s.absolute_model_check_pass = res.absPass;
s.quat_norm_error_collocation = res.qNormColl;
s.quat_norm_error_propagation = res.qNormProp;
s.burns = rmfield(res.burns, {'n0', 'n1'});
s.thruster_impulse_Ns = res.thrusterImpulse';
s.reference_fit_error_m = res.refFitErr;
s.sun = struct('anti_sun_angle_min_deg', min(res.sun.antiDeg), ...
    'anti_sun_angle_max_deg', max(res.sun.antiDeg), ...
    'max_violation_deg', res.sun.maxViolDeg, 'pass', res.sun.pass, ...
    'time_outside_s', res.sun.timeOutside, 'outside_intervals', {res.sun.intervals}, ...
    'collocation_max_violation_deg', res.sun.collMaxViolDeg, ...
    'formation_centre_error_rad', res.sun.centreErr, ...
    'active_bound', res.sun.activeBound, 'active_points', {res.sun.active});
if isfield(res, 'multistart') && ~isempty(res.multistart)
    s.multistart = res.multistart;
    writeText(fullfile(outDir, 'multistart.json'), jsonencode(res.multistart, 'PrettyPrint', true));
end
writeText(fullfile(outDir, 'summary.json'), jsonencode(s, 'PrettyPrint', true));

settings = rmfield(p, {'thr'});
settings.thruster_pos = p.thr.pos; settings.thruster_dir = p.thr.dir;
settings.thruster_label = p.thr.label;
writeText(fullfile(outDir, 'settings.json'), jsonencode(settings, 'PrettyPrint', true));

sol = res.sol;
lbl = cellstr("T" + (1:p.nThr) + "_N");
writetable(array2table([sol.tNode(1:end-1)', sol.tNode(2:end)', sol.T'], ...
    'VariableNames', [{'t_start_s', 't_end_s'}, lbl]), fullfile(outDir, 'thrust_history.csv'));
xn = {'qw', 'qx', 'qy', 'qz', 'wx', 'wy', 'wz', 'rx', 'ry', 'rz', 'vx', 'vy', 'vz'};
writetable(array2table([res.prop.t', res.prop.x', res.ref', rad2deg(res.sun.theta)'], 'VariableNames', ...
    [{'t_s'}, xn, {'rbf_x', 'rbf_y', 'rbf_z', 'vbf_x', 'vbf_y', 'vbf_z', 'sun_angle_deg'}]), ...
    fullfile(outDir, 'propagated_states.csv'));
writetable(array2table([sol.tNode', sol.X'], 'VariableNames', [{'t_s'}, xn]), ...
    fullfile(outDir, 'collocation_states.csv'));
end

function writeText(file, txt)
fid = fopen(file, 'w'); fprintf(fid, '%s', txt); fclose(fid);
end
