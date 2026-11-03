%% Problem 3 phase sweep, Sec. 5.2 validation (run after main.m)
%  Problem 3 with the arrival phase phiTau fixed at 16 equally spaced
%  values, each from the common initial guess. The multi-start optimum of
%  main.m should match the sweep minimum (J3 <= J2 must hold). The
%  prescribed arrival phase of Problem 2 is solved as well: Problem 3
%  with phiTau = phi0 + rate*tau is Problem 2 (reduction check).
projectRoot = fileparts(mfilename('fullpath'));
addpath(projectRoot);
addpath(genpath(fullfile(projectRoot, 'src')));
setupCasadi();

p = lifeParameters();
p.ocp.ipopt.print_level = 0;
ref = propagateReference(p, p.tf);
ocp = buildOcp(p, ref, 3);
phiPresc = p.target.phase0 + p.target.rate*p.tf;    % Problem 2 at tau = tf
phis = [-pi + (0:15)*pi/8, phiPresc];
sweep = struct('phiTauDeg', {}, 'J_Ns', {}, 'tau_s', {}, 'status', {}, ...
               'iterations', {}, 'time_s', {}, 'prescribed', {});
for k = 1:numel(phis)
    o = ocp;
    o.lbx(o.ix.phi) = phis(k)/pi; o.ubx(o.ix.phi) = phis(k)/pi;
    [s, in] = solveOcp(o, initialGuess(o, p, phis(k)), p);
    sweep(k) = struct('phiTauDeg', rad2deg(phis(k)), 'J_Ns', s.J, 'tau_s', s.tau, ...
        'status', in.status, 'iterations', in.iter, 'time_s', in.wallTime, ...
        'prescribed', k == numel(phis));
    fprintf('phiTau %7.2f deg: J %.4f N s, tau %.1f s, %s (%d it., %.1f s)\n', ...
            sweep(k).phiTauDeg, s.J, s.tau, in.status, in.iter, in.wallTime);
end

% Checks against the saved Problems 2 and 3
d3 = load(fullfile(projectRoot, 'results', 'problem3', 'solution.mat'), 'res');
d2 = load(fullfile(projectRoot, 'results', 'problem2', 'solution.mat'), 'res');
ok = strcmp({sweep.status}, 'Solve_Succeeded') & ~[sweep.prescribed];
[Jmin, iMin] = min([sweep(ok).J_Ns]);
okIdx = find(ok);
chk.J3_Ns = d3.res.J;
chk.phiTau3_deg = rad2deg(d3.res.sol.phiTau);
chk.sweep_min_J_Ns = Jmin;
chk.sweep_min_phiTau_deg = sweep(okIdx(iMin)).phiTauDeg;
chk.J3_le_sweep_min = d3.res.J <= Jmin*(1 + 1e-6);
chk.J2_Ns = d2.res.J;
chk.prescribed_phase_J_Ns = sweep(end).J_Ns;
chk.reduction_difference_Ns = sweep(end).J_Ns - d2.res.J;
fprintf(['\nPhase sweep: min J %.4f N s at phiTau %.1f deg; multi-start J3 %.4f N s at %.2f deg ' ...
         '-> J3 <= sweep minimum: %d\n'], Jmin, chk.sweep_min_phiTau_deg, chk.J3_Ns, ...
        chk.phiTau3_deg, chk.J3_le_sweep_min);
fprintf('Reduction check (phiTau prescribed, %.1f deg): J %.4f N s vs Problem 2 J2 %.4f N s (diff %+.4f)\n', ...
        rad2deg(phiPresc), sweep(end).J_Ns, d2.res.J, chk.reduction_difference_Ns);

outDir = fullfile(projectRoot, 'results', 'problem3');
fid = fopen(fullfile(outDir, 'phase_sweep.json'), 'w');
fprintf(fid, '%s', jsonencode(struct('sweep', sweep, 'checks', chk), 'PrettyPrint', true));
fclose(fid);
fig = plotPhaseSweep(sweep, d3.res, d2.res, p);
exportgraphics(fig, fullfile(outDir, 'phase_sweep.png'), 'Resolution', 200);
cmpDir = fullfile(projectRoot, 'results', 'comparison');
if isfolder(cmpDir)
    exportgraphics(fig, fullfile(cmpDir, 'phase_sweep.png'), 'Resolution', 200);
end
