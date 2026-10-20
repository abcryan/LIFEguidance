function [levels, final] = runMeshRefinement(p, ref, logDir)
%RUNMESHREFINEMENT Solve on the initial uniform mesh, double it once (the
%   "doubling N" check of Sec. 4.3), then refine nonuniformly near thrust
%   switches and where the local discretization error is large.
%   A level is accepted when
%     - IPOPT returns Solve_Succeeded,
%     - the independent propagation passes (59)-(60),
%     - total impulse changes < p.mesh.relTolJ from the previous level.
%   Refinement continues while some burn's inertial force/torque integral,
%   duration or peak thrust still changes (near-impulsive burns shorten at
%   every bisection), up to the caps. It stops early, keeping the last
%   accepted level, if a warm-started level does not reach Solve_Succeeded
%   or jumps to another local optimum (|dJ| above the tolerance).
%   final: index into levels of the reported solution.
sigma = ones(1, p.ocp.N0)/p.ocp.N0;
prev = [];
levels = struct([]);
final = 0;
for L = 0:p.mesh.maxLevels
    warm = ~isempty(prev);
    ocp = buildOcp(p, ref, sigma, warm);
    if warm, w0 = warmStart(ocp, prev, p); else, w0 = initialGuess(ocp, p); end
    logFile = fullfile(logDir, sprintf('ipopt_level%d.log', L));
    [sol, info] = solveOcp(ocp, w0, p, logFile);
    info.retried = false;
    if warm && ~info.success                % continue once from where it stopped
        [sol, info] = solveOcp(ocp, sol.w, p, strrep(logFile, '.log', '_retry.log'));
        info.retried = true;
    end
    prop = propagateRelative(sol, p, 9);
    term = terminalErrors(prop.x(:,end), p);
    loc = localErrors(sol, ref, p);
    burns = detectBurns(sol, prop, p);

    lv.level = L; lv.N = sol.N; lv.sol = sol; lv.info = info;
    lv.term = term; lv.loc = loc; lv.burns = burns;
    lv.uniform = all(abs(sol.sigma - 1/sol.N) < 1e-12);
    if isempty(prev)
        lv.dJ = NaN; lv.burnCmp = [];
    else
        lv.dJ = abs(sol.J - prev.J)/sol.J;
        lv.burnCmp = compareBurns(burns, levels(end).burns, p);
    end
    lv.accepted = info.success && term.passAll && lv.dJ < p.mesh.relTolJ;
    lv.qNormColl = max(abs(vecnorm([sol.X(1:4,:), sol.XC(1:4,:)]) - 1));
    lv.qNormProp = max(abs(vecnorm(prop.x(1:4,:)) - 1));
    levels = [levels, lv]; %#ok<AGROW>

    fprintf(['level %d  N %4d  %-16s iter %4d  %6.1f s | tau %8.2f s  ' ...
             'J %9.4f N s  dJ %8.2e | r %.1e v %.1e th %.1e w %.1e  accepted %d\n'], ...
            L, sol.N, info.status, info.iter, info.wallTime, sol.tau, sol.J, ...
            lv.dJ, term.r, term.v, term.thetaAcos, term.w, lv.accepted);

    if L == 0
        sigma = repelem(sol.sigma/2, 2);    % doubling check
        prev = sol;
        continue
    end
    if ~lv.accepted
        fprintf('level %d not accepted: refinement stopped\n', L);
        break
    end
    final = numel(levels);
    if lv.burnCmp.ok
        fprintf('burn features converged at level %d\n', L);
        break
    end
    [sigma, marked] = refineMesh(sol, loc, p);
    if ~any(marked)
        fprintf('no interval marked for refinement (or N cap reached)\n');
        break
    end
    prev = sol;
end
if final == 0
    warning('runMeshRefinement:noAccepted', 'No level was accepted; reporting the last one.');
    final = numel(levels);
end
end
