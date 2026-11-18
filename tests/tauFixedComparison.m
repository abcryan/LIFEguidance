%% Problem 3 with the arrival time tau free vs fixed to tf (benchmark)
%  Collector 1, the 8-start multi-start of main.m with tau free and with
%  tau fixed to tf (lbx = ubx = 1, which IPOPT removes from the problem),
%  on the same compiled solver; run order free, fixed, free (timing noise).
%  Prints iterations, time per iteration and the best cost of each run, and
%  writes results/problem3_tau_fixed/tauFixedComparison.mat. ~12 min.
%  Not a unit test (the name does not end in Test, so runtests skips it).
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root); addpath(genpath(fullfile(root, 'src'))); setupCasadi();
p = lifeParameters();
p.ocp.ipopt.print_timing_statistics = 'yes';
ref = propagateReference(p, p.tf);
tB = tic; ocp = buildOcp(p, ref, 3); fprintf('build %.1f s\n', toc(tB));
ocpF = ocp; ocpF.lbx(ocp.ix.tau) = 1; ocpF.ubx(ocp.ix.tau) = 1;
runs = {'free', ocp; 'fixed', ocpF; 'free (repeat)', ocp};
for r = 1:size(runs, 1)
    fprintf('\n######## RUN %s ########\n', runs{r, 1});
    tR = tic;
    [sol, info, st] = solveMultiStart(runs{r, 2}, p);
    R(r) = struct('name', runs{r, 1}, 'starts', st, 'J', sol.J, 'phiTauDeg', rad2deg(sol.phiTau), ...
        'tau', sol.tau, 'tauMultiplier', info.tauMultiplier, 'wall', toc(tR)); %#ok<SAGROW>
end
outDir = fullfile(root, 'results', 'problem3_tau_fixed');
if ~isfolder(outDir), mkdir(outDir); end
save(fullfile(outDir, 'tauFixedComparison.mat'), 'R');
fprintf('\n######## SUMMARY ########\n');
for r = 1:numel(R)
    st = R(r).starts; ok = strcmp({st.status}, 'Solve_Succeeded');
    it = [st.iterations]; tm = [st.time_s];
    fprintf('%-14s best J %.4f at %.2f deg, tau %.1f | total %.1f s, %d it., %.3f s/it | succeeded %d/8 | per start (it/s): %s\n', ...
        R(r).name, R(r).J, R(r).phiTauDeg, R(r).tau, sum(tm), sum(it), sum(tm)/sum(it), sum(ok), ...
        strjoin(arrayfun(@(k) sprintf('%d/%.1f', it(k), tm(k)), 1:numel(it), 'UniformOutput', false), ' '));
end
