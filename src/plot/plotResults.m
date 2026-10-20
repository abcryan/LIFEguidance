function figs = plotResults(res, p, outDir)
%PLOTRESULTS Static figures of Sec. 4.4, saved as PNG to outDir and left
%   open (physical time throughout).
list = {
    'tracking_errors',     @() plotTracking(res, p)
    'control_wrench',      @() plotWrench(res, p)
    'thrust_histories',    @() plotThrusters(res, p)
    'thrust_burn_zoom',    @() plotThrusterBurns(res, p)
    'thrust_heatmap',      @() plotHeatmap(res, p)
    'impulse',             @() plotImpulse(res, p)
    'propagation_discrepancy', @() plotDiscrepancy(res, p)
    'trajectory_formation', @() plotTrajectoryFormation(res, p)
    'trajectory_SER',      @() plotTrajectorySER(res, p)
    'halo_orbit',          @() plotHaloOrbit(res, p)
    };
for k = 1:size(list, 1)
    figs.(list{k, 1}) = list{k, 2}();
    exportgraphics(figs.(list{k, 1}), fullfile(outDir, [list{k, 1} '.png']), 'Resolution', 200);
end
end
