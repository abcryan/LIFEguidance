function fig = plotPhaseSweep(sweep, res3, res2, ~)
%PLOTPHASESWEEP Sec. 5.3, Problem 3: J_imp against the fixed arrival phase
%   (phase sweep), the multi-start results and the optimum phiTau*, and the
%   prescribed arrival phase of Problem 2 with its cost J2.
s = plotStyle();
fig = newFigure('phaseSweep', [50 50 1100 550]);
ax = axes(fig); hold(ax, 'on'); grid(ax, 'on');
sw = sweep(~[sweep.prescribed]);
ok = strcmp({sw.status}, 'Solve_Succeeded');
x = [sw(ok).phiTauDeg]; y = [sw(ok).J_Ns];
if any(x == -180), x = [x, 180]; y = [y, y(x == -180)]; end   % 2*pi periodic
plot(ax, x, y, 'o-', 'Color', [0.4 0.4 0.4], ...
     'MarkerFaceColor', [0.6 0.6 0.6], 'DisplayName', 'sweep (\phi_\tau fixed)');
if any(~ok)
    plot(ax, [sw(~ok).phiTauDeg], [sw(~ok).J_Ns], 'rx', 'MarkerSize', 10, 'LineWidth', 2, ...
         'DisplayName', 'sweep, not converged');
end
ms = res3.multistart;
conv = strcmp({ms.status}, 'Solve_Succeeded');
plot(ax, [ms(conv).phiTauDeg], [ms(conv).J_Ns], 'd', 'MarkerSize', 8, 'Color', s.prob(3,:), ...
     'MarkerFaceColor', 'w', 'LineWidth', 1.5, 'DisplayName', ...
     sprintf('multi-start results (%d of %d converged)', nnz(conv), numel(ms)));
plot(ax, rad2deg(res3.sol.phiTau), res3.J, 'p', 'MarkerSize', 16, 'Color', 'k', ...
     'MarkerFaceColor', s.prob(3,:), 'DisplayName', sprintf('Problem 3: \\phi_\\tau^* = %.1f\\circ, J_3 = %.2f N s', ...
     rad2deg(res3.sol.phiTau), res3.J));
phiP = rad2deg(res2.sol.phiTau);
xline(ax, phiP, '--', 'Color', s.prob(2,:), 'LineWidth', 1.5, 'HandleVisibility', 'off');
plot(ax, phiP, res2.J, 's', 'MarkerSize', 10, 'Color', 'k', 'MarkerFaceColor', s.prob(2,:), ...
     'DisplayName', sprintf('Problem 2: prescribed \\phi_\\tau = %.0f\\circ, J_2 = %.2f N s', phiP, res2.J));
if any([sweep.prescribed])
    pr = sweep([sweep.prescribed]);
    plot(ax, pr.phiTauDeg, pr.J_Ns, '+', 'MarkerSize', 12, 'LineWidth', 2, 'Color', s.prob(2,:), ...
         'DisplayName', sprintf('Problem 3 at the prescribed phase: %.2f N s', pr.J_Ns));
end
xlim(ax, [-180 180]); xticks(ax, -180:45:180);
xlabel(ax, 'Arrival phase \phi_\tau (deg)'); ylabel(ax, 'Total impulse J_{imp} (N s)');
legend(ax, 'Location', 'north');
title(ax, sprintf('Problem 3 phase sweep (%d fixed phases), multi-start optimum and Problem 2', ...
      numel(sw)), 'FontWeight', 'normal');
end
