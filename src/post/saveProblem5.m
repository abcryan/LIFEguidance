function saveProblem5(F, M, ref, p, outDir)
%SAVEPROBLEM5 Machine-readable outputs of Problem 5 in outDir:
%   solution.mat (fleet: per-spacecraft results, starts, sweep, collision
%   multipliers; metrics; p), summary.json, starts.json, dmin_sweep.json,
%   and scI/ per spacecraft (saveSpacecraft).
if ~isfolder(outDir), mkdir(outDir); end
fleet = F; metrics = M;
save(fullfile(outDir, 'solution.mat'), 'fleet', 'metrics', 'p', '-v7.3');
s.problem = 5;
s.status = F.info.status; s.iterations = F.info.iter; s.solve_time_s = F.info.wallTime;
s.nlp_build_time_s = F.buildTime; s.validation_time_s = F.validationTime;
s.N = p.ocp.N; s.dmin_m = p.collision.dmin; s.nlp_margin_m = p.collision.margin;
s.tau_s = M.tau; s.phi_tau_deg = rad2deg(M.phiTau);
s.total_impulse_Ns = M.J5; s.propellant_kg = M.dm;
s.J4_Ns = M.J4; s.Pi5_Ns = M.Pi5; s.Pi5_relative = M.Pi5rel;
s.delta_tau_s = M.dTau; s.delta_phi_tau_deg = M.dPhiDeg;
s.dJ5_ddmin_multipliers_Ns_per_m = M.dJddmin;
s.dJ5_ddmin_central_differences_Ns_per_m = M.dJddminFD;
s.dJ5_ddmin_sweep_slope_Ns_per_m = M.sweepSlopeAtDmin;
s.crossing_pairs = F.crossing;
s.active_pairs = M.activePairs;
s.active_intervals = M.activeIntervals;
s.collocation = M.colloc;
s.min_separation = struct('m', M.minSep, 'pair', M.minSepPair, 't_s', M.minSepTime);
s.separation_minima = M.minima;
s.formations = M.feasibility;
for i = numel(F.res):-1:1
    r = F.res{i};
    sc(i) = struct('spacecraft', i, 'name', p.fleet.name(i), 'J_Ns', M.Ji(i), 'J4_Ns', M.J4i(i), ...
        'propellant_kg', M.dmi(i), 'max_thrust_N', r.maxThrust, ...
        'sun_max_violation_deg', r.sun.maxViolDeg, 'terminal', rmfield(r.term, {'pass', 'passAll'}), ...
        'terminal_pass', r.term.passAll, 'absolute_model_check_pass', r.absPass, 'accepted', r.accepted);
end
s.spacecraft = sc;
s.checks = M.checks; s.check_tolerance_J_Ns = M.tolJ;
s.expectations = M.expectations;
s.starts = F.starts; s.dmin_sweep = F.sweep; s.finite_differences = F.fd;
s.low_dmin_from_problem4 = F.lowFromP4;
writeText(fullfile(outDir, 'summary.json'), jsonencode(s, 'PrettyPrint', true));
writeText(fullfile(outDir, 'starts.json'), jsonencode(F.starts, 'PrettyPrint', true));
writeText(fullfile(outDir, 'dmin_sweep.json'), jsonencode(F.sweep, 'PrettyPrint', true));
saveSpacecraft(F.res, ref, outDir);
end

function writeText(file, txt)
fid = fopen(file, 'w'); fprintf(fid, '%s', txt); fclose(fid);
end
