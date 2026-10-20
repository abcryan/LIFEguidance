function saveResults(res, levels, ref, p, outDir)
%SAVERESULTS Machine-readable outputs for reproducibility (Sec. 4.4).
%   solution.mat (v7.3/HDF5), summary.json, settings.json and CSV files.
ref = rmfield(ref, 'eval');             % function handles are not portable
save(fullfile(outDir, 'solution.mat'), 'res', 'levels', 'ref', 'p', '-v7.3');

s.status = res.info.status;
s.iterations = res.info.iter;
s.accepted = res.accepted;
s.N = res.sol.N;
s.tau_s = res.sol.tau;
s.tau_at_upper_bound = res.info.tauAtUpper;
s.tau_bound_multiplier = res.info.tauMultiplier;
s.total_impulse_Ns = res.J;
s.propellant_kg = res.dm;
s.max_thrust_N = res.maxThrust;
s.thrust_bound_violation_N = res.boundViolation;
s.terminal = rmfield(res.term, {'pass', 'passAll'});
s.terminal_pass = res.term.pass;
s.terminal_pass_all = res.term.passAll;
s.doubling_dJ = res.doublingDJ;
s.last_level_dJ = res.history.dJ(res.final);
s.prop_vs_collocation = res.disc.max;
s.absolute_model_check = res.absDisc.max;
s.absolute_model_check_pass = res.absPass;
s.quat_norm_error_collocation = res.qNormColl;
s.quat_norm_error_propagation = res.qNormProp;
s.burns_converged = res.burnsConverged;
s.burns = rmfield(res.burns, {'n0', 'n1'});
s.thruster_impulse_Ns = res.thrusterImpulse';
s.reference_fit_error_m = res.refFitErr;
s.levels = table2struct(res.history);
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
writetable(array2table([res.prop.t', res.prop.x', res.ref'], 'VariableNames', ...
    [{'t_s'}, xn, {'rbf_x', 'rbf_y', 'rbf_z', 'vbf_x', 'vbf_y', 'vbf_z'}]), ...
    fullfile(outDir, 'propagated_states.csv'));
writetable(array2table([sol.tNode', sol.X'], 'VariableNames', [{'t_s'}, xn]), ...
    fullfile(outDir, 'collocation_states.csv'));
writetable(res.history, fullfile(outDir, 'mesh_convergence.csv'));
end

function writeText(file, txt)
fid = fopen(file, 'w'); fprintf(fid, '%s', txt); fclose(fid);
end
