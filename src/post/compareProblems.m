function cmp = compareProblems(resultsDir, p)
%COMPAREPROBLEMS Comparison of the saved solutions (Sec. 5.3): table of
%   tau*, J_imp, dm, phi_tau*, active constraints, Sun-angle violation and
%   solve time, the consistency checks of Sec. 4.5, and the comparison
%   figures. Reads results/problemK/solution.mat for every K present (each
%   with its own parameters, R{i}.p, which carry the solution's formation
%   phase) and writes results/comparison/. The Problem 3 phase sweep
%   (phaseSweep.m) is included when present. Problem 4 enters the table
%   and the figures through collector 1, the spacecraft of Problems 1-3
%   (at the shared tau* and phiTau*), and adds a fleet block: the total
%   and per-spacecraft impulse, Pi4 and the fleet checks.
R = {};
for k = 1:3
    f = fullfile(resultsDir, sprintf('problem%d', k), 'solution.mat');
    if isfile(f)
        d = load(f, 'res', 'p');
        d.res.p = d.p;
        R{end + 1} = d.res; %#ok<AGROW>
    end
end
f4 = fullfile(resultsDir, 'problem4', 'solution.mat');
fleet = [];
if isfile(f4)
    d = load(f4, 'fleet', 'metrics', 'sweep');
    R{end + 1} = d.fleet.res{1};        % collector 1, with its own p
    fleet = d;
end
if numel(R) < 2
    fprintf('Comparison skipped: fewer than two problems in %s\n', resultsDir);
    cmp = [];
    return
end
outDir = fullfile(resultsDir, 'comparison');
if ~isfolder(outDir), mkdir(outDir); end

for i = 1:numel(R)
    r = R{i};
    c(i).problem = r.problem; %#ok<AGROW>
    c(i).status = r.info.status;
    c(i).iterations = r.info.iter;
    c(i).solve_time_s = r.info.wallTime;
    c(i).total_solve_time_s = r.info.wallTime;  % Problem 3: all multi-start solves
    if isfield(r.info, 'multiStartTime'), c(i).total_solve_time_s = r.info.multiStartTime; end
    if r.problem == 4                   % joint starts and the sweep
        c(i).total_solve_time_s = fleet.fleet.jointTime + fleet.sweep.wallTime;
    end
    c(i).tau_s = r.sol.tau;
    c(i).tau_at_upper_bound = r.info.tauAtUpper;
    c(i).J_imp_Ns = r.J;
    c(i).propellant_kg = r.dm;
    c(i).phi_tau_deg = rad2deg(angle(exp(1i*r.sol.phiTau)));
    c(i).phi_tau_free = r.problem >= 3;
    c(i).max_thrust_N = r.maxThrust;
    c(i).thrust_bound_active = r.maxThrust >= (1 - 1e-6)*p.Tmax;
    c(i).sun_max_violation_deg = r.sun.maxViolDeg;
    c(i).sun_time_outside_s = r.sun.timeOutside;
    c(i).anti_sun_min_deg = min(r.sun.antiDeg);
    c(i).anti_sun_max_deg = max(r.sun.antiDeg);
    c(i).sun_active_points = numel(r.sun.active);
    c(i).sun_active_times_s = strjoin(compose('%.0f', [r.sun.active.t]), ' ');
    c(i).sun_active_bound = r.sun.activeBound;
    c(i).arrival_checks_pass = r.term.passAll;
    c(i).accepted = r.accepted;
end
cmp.table = struct2table(c);

% Consistency checks (Sec. 4.5); local optima, so these are checks only
P = [c.problem];
J = @(k) c(P == k).J_imp_Ns;
cmp.checks = struct();
if all(ismember([1 2], P))
    cmp.checks.nesting_J1_le_J2 = J(1) <= J(2)*(1 + 1e-9);
    cmp.checks.problem1_violates_sun = c(P == 1).sun_max_violation_deg > p.gate.sunDeg;
    cmp.checks.problem2_satisfies_sun = c(P == 2).sun_max_violation_deg <= p.gate.sunDeg;
    cmp.checks.J2_minus_J1_Ns = J(2) - J(1);
end
if all(ismember([2 3], P))
    cmp.checks.nesting_J3_le_J2 = J(3) <= J(2)*(1 + 1e-9);
    cmp.checks.problem3_satisfies_sun = c(P == 3).sun_max_violation_deg <= p.gate.sunDeg;
    cmp.checks.J3_minus_J2_Ns = J(3) - J(2);
end
sweepFile = fullfile(resultsDir, 'problem3', 'phase_sweep.json');
if any(P == 3) && isfile(sweepFile)
    sw = jsondecode(fileread(sweepFile));
    cmp.checks.J3_le_phase_sweep_min = J(3) <= sw.checks.sweep_min_J_Ns*(1 + 1e-6);
    cmp.checks.phase_sweep_min_Ns = sw.checks.sweep_min_J_Ns;
    cmp.checks.prescribed_phase_minus_J2_Ns = sw.checks.reduction_difference_Ns;
    png = fullfile(resultsDir, 'problem3', 'phase_sweep.png');
    if isfile(png), copyfile(png, fullfile(outDir, 'phase_sweep.png')); end
end

cmp.fleet = [];
if any(P == 4)
    M = fleet.metrics;
    cmp.checks.problem4_collector1_satisfies_sun = c(P == 4).sun_max_violation_deg <= p.gate.sunDeg;
    if any(P == 3)
        cmp.checks.collector1_J4_minus_J3_Ns = J(4) - J(3);
        cmp.checks.problem3_J3_minus_fleet_J3_1_Ns = J(3) - M.J3i(1);
    end
    fn = fieldnames(M.checks);
    for k = 1:numel(fn)
        cmp.checks.(['fleet_' fn{k}]) = M.checks.(fn{k});
    end
    cmp.fleet = struct('J4_Ns', M.J4, 'propellant_kg', M.dm, 'tau_s', M.tau, ...
        'phi_tau_deg', rad2deg(M.phiTau), 'J_per_spacecraft_Ns', M.Ji, 'share', M.share, ...
        'J3_per_spacecraft_Ns', M.J3i, 'Pi4_Ns', M.Pi4, 'Pi4_relative', M.Pi4rel, ...
        'combiner_share', M.share(5), 'min_distance_m', M.dist.min, ...
        'sweep_time_s', fleet.sweep.wallTime, 'joint_time_s', fleet.fleet.jointTime);
end

writetable(cmp.table, fullfile(outDir, 'comparison.csv'));
fid = fopen(fullfile(outDir, 'comparison.json'), 'w');
fprintf(fid, '%s', jsonencode(struct('problems', c, 'fleet', cmp.fleet, 'checks', cmp.checks), ...
    'PrettyPrint', true));
fclose(fid);
printComparison(c, cmp.checks, cmp.fleet);

cmp.figs = plotComparison(R, p, outDir);
end

function printComparison(c, checks, fleet)
fprintf('\n===== Comparison of the guidance problems (Problem 4: collector 1 of the fleet) =====\n');
row = @(name, fmt, v) fprintf(['%-28s' repmat(sprintf('  %%%ss', fmt), 1, numel(v)) '\n'], name, v{:});
str = @(f) arrayfun(f, c, 'UniformOutput', false);
row('Problem', '16', str(@(x) sprintf('%d', x.problem)));
row('Status', '16', str(@(x) x.status(1:min(end, 16))));
row('Iterations / solve time', '16', str(@(x) sprintf('%d / %.1f s', x.iterations, x.solve_time_s)));
row('Total solve time (s)', '16', str(@(x) sprintf('%.1f', x.total_solve_time_s)));
row('Arrival time tau* (s)', '16', str(@(x) sprintf('%.1f%s', x.tau_s, star(x.tau_at_upper_bound))));
row('Total impulse J (N s)', '16', str(@(x) sprintf('%.3f', x.J_imp_Ns)));
row('Propellant dm (kg)', '16', str(@(x) sprintf('%.4f', x.propellant_kg)));
row('Arrival phase phi_tau (deg)', '16', str(@(x) sprintf('%.2f%s', x.phi_tau_deg, phaseTag(x.phi_tau_free))));
row('Max thrust (N)', '16', str(@(x) sprintf('%.3f%s', x.max_thrust_N, star(x.thrust_bound_active))));
row('Anti-Sun angle range (deg)', '16', str(@(x) sprintf('[%.1f, %.1f]', x.anti_sun_min_deg, x.anti_sun_max_deg)));
row('Max Sun violation (deg)', '16', str(@(x) sprintf('%.3f', x.sun_max_violation_deg)));
row('Time outside > gate (s)', '16', str(@(x) sprintf('%.0f', x.sun_time_outside_s)));
row('Sun constraint active (pts)', '16', str(@(x) sprintf('%d', x.sun_active_points)));
row('Arrival checks / accepted', '16', str(@(x) sprintf('%d / %d', x.arrival_checks_pass, x.accepted)));
fprintf('(* bound active; p prescribed, f free)\n');
if ~isempty(fleet)
    fprintf(['Problem 4 fleet: J4* %.3f N s (dm %.4f kg), per spacecraft %s N s, combiner %.1f%%, ' ...
             'sum J*_3,i %.3f N s, Pi4 %.3f N s (%.2f%%), min. distance %.1f m\n'], fleet.J4_Ns, ...
            fleet.propellant_kg, mat2str(fleet.J_per_spacecraft_Ns, 6), 100*fleet.combiner_share, ...
            sum(fleet.J3_per_spacecraft_Ns), fleet.Pi4_Ns, 100*fleet.Pi4_relative, fleet.min_distance_m);
end
fn = fieldnames(checks);
for k = 1:numel(fn)
    fprintf('Check %-26s: %s\n', fn{k}, num2str(checks.(fn{k})));
end
fprintf('\n');
end

function c = star(b)
if b, c = '*'; else, c = ''; end
end

function c = phaseTag(free)
if free, c = ' f'; else, c = ' p'; end
end
