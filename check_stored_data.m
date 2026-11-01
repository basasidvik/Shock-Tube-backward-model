% CHECK_STORED_DATA  Show what the generated files actually contain.
%
% Run after step1 (and step2). Lists every full .mat file (ALL grid points)
% and every sensor CSV (only the sensor columns). Reads the settings from
% shock_tube_config.m.

clc; clear;
cfg = shock_tube_config();

mats = dir([cfg.filePrefix '_t*.mat']);
fprintf('FULL datasets (.mat) - all grid points are stored here\n');
fprintf('%-34s %7s %17s %12s %8s %7s\n', 'file', 't', 'X (cases x cols)', 'Y', 'MB', 'used?');
for k = 1:numel(mats)
    S = load(mats(k).name);  d = S.data;
    if any(abs(cfg.useTimes - d.opts.tEval) < 1e-9), used = 'yes'; else, used = 'no'; end
    fprintf('%-34s %7.3f %9d x %-6d %6d x %d %8.0f %7s\n', mats(k).name, d.opts.tEval, ...
            size(d.X, 1), size(d.X, 2), size(d.Y, 1), size(d.Y, 2), ...
            mats(k).bytes / 1e6, used);
end
if ~isempty(mats)
    fprintf('-> each X row = %d grid points x 3 variables (rho | u | p) = %d numbers\n', ...
            numel(d.x), size(d.X, 2));
else
    fprintf('(none found - run step1_generate_datasets first)\n');
end

csvs = dir([cfg.sensorPrefix '_*.csv']);
fprintf('\nSENSOR files (.csv) - only the sensor positions are stored here\n');
fprintf('%-44s %8s %10s\n', 'file', 'columns', 'MB');
for k = 1:numel(csvs)
    fid = fopen(csvs(k).name);  h = fgetl(fid);  fclose(fid);
    fprintf('%-44s %8d %10.2f\n', csvs(k).name, numel(strfind(h, ',')) + 1, csvs(k).bytes / 1e6);
end
nS = numel(cfg.sensorX);  nV = numel(cfg.sensorVars);  nT = numel(cfg.useTimes);
fprintf('-> a single-time X file has 1 case_id + %d variables x %d sensors = %d columns\n', ...
        nV, nS, 1 + nV * nS);
fprintf('   the alltimes file has 1 + %d x %d x %d times = %d columns\n', ...
        nV, nS, nT, 1 + nV * nS * nT);
