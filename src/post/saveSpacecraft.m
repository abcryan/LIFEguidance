function saveSpacecraft(res, ref, outDir)
%SAVESPACECRAFT Single-spacecraft outputs of a fleet solution (Problems 4
%   and 5): for every spacecraft outDir/scI/ with saveResults and the
%   figures of plotResults, written in parallel (figures closed).
nS = numel(res);
parfor i = 1:nS
    d = fullfile(outDir, sprintf('sc%d', i));
    saveResults(res{i}, ref, res{i}.p, d);
    figs = plotResults(res{i}, res{i}.p, d);
    structfun(@close, figs);
end
end
