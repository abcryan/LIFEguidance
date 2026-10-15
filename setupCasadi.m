function setupCasadi()
%SETUPCASADI Add CasADi to the path, downloading it into external/ if missing.
version = '3.8.1';
target = fullfile(fileparts(mfilename('fullpath')), 'external', 'casadi');

if ~isfolder(fullfile(target, '+casadi'))
    assets = containers.Map({'maca64', 'maci64', 'glnxa64', 'win64'}, ...
        {'osx-arm64', 'osx-x86_64', 'linux64', 'windows64'});
    arch = computer('arch');
    if ~isKey(assets, arch)
        error('setupCasadi:unsupported', 'No CasADi build for %s.', arch);
    end
    url = sprintf(['https://github.com/casadi/casadi/releases/download/' ...
                   '%s/casadi-%s-%s-matlab2018b.zip'], version, version, assets(arch));
    fprintf('Downloading CasADi %s (~100 MB) ...\n', version);
    zipFile = [tempname '.zip'];
    websave(zipFile, url);
    unzip(zipFile, target);
    delete(zipFile);
end
addpath(target);
end
