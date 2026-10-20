function plotResults(res, levels, p, outDir)
%PLOTRESULTS Static figures of Sec. 4.4 (PNG), physical time throughout.
figs = {
    'tracking_errors',     @() plotTracking(res, p)
    'control_wrench',      @() plotWrench(res, p)
    'thrust_histories',    @() plotThrusters(res, p)
    'thrust_burn_zoom',    @() plotThrusterBurns(res, p)
    'thrust_heatmap',      @() plotHeatmap(res, p)
    'impulse',             @() plotImpulse(res, p)
    'propagation_discrepancy', @() plotDiscrepancy(res, p)
    'mesh_convergence',    @() plotMeshConvergence(levels, p)
    'trajectory_formation', @() plotTrajectoryFormation(res, p)
    'trajectory_SER',      @() plotTrajectorySER(res, p)
    };
for k = 1:size(figs, 1)
    fig = figs{k, 2}();
    exportgraphics(fig, fullfile(outDir, [figs{k, 1} '.png']), 'Resolution', 200);
    close(fig);
end
end
