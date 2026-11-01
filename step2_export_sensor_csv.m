% STEP2_EXPORT_SENSOR_CSV  Use only the chosen times (cfg.useTimes) and write
% labelled CSVs containing just the sensor positions (cfg.sensorX).
%
% The full .mat files from step1 are only read, never changed, so you can
% rerun this step with other sensor positions, noise, or times without
% regenerating the datasets. Settings: shock_tube_config.m

clc; clear;
cfg = shock_tube_config();

% ---- pick the files for the times that are actually used ------------------
files = cell(numel(cfg.useTimes), 1);
for i = 1:numel(cfg.useTimes)
    files{i} = sprintf('%s_t%.3f.mat', cfg.filePrefix, cfg.useTimes(i));
    if exist(files{i}, 'file') ~= 2
        error('%s not found. Run step1_generate_datasets first, and make sure %.3f is in cfg.allTimes.', ...
              files{i}, cfg.useTimes(i));
    end
end

% ---- sensor CSVs (the training data for the ML model) ---------------------
sopts.mode     = cfg.sensorMode;
sopts.vars     = cfg.sensorVars;
sopts.noiseStd = cfg.sensorNoise;
sopts.prefix   = cfg.sensorPrefix;
export_sensor_csv(files, cfg.sensorX, sopts);

% ---- optional: all grid points to labelled CSV (very large) ---------------
if cfg.exportFullCsv
    for it = 1:numel(files)
        S = load(files{it});  data = S.data;
        nx     = numel(data.x);
        caseID = (1:size(data.X,1))';
        xNames = [compose('rho_%03d',1:nx), compose('u_%03d',1:nx), compose('p_%03d',1:nx)];
        Xtbl   = array2table(double(data.X), 'VariableNames', xNames);
        Xtbl   = addvars(Xtbl, caseID, 'Before', 1, 'NewVariableNames', 'case_id');
        Ytbl   = array2table(data.Y, 'VariableNames', data.targetNames);
        Ytbl   = addvars(Ytbl, caseID, 'Before', 1, 'NewVariableNames', 'case_id');
        base = files{it}(1:end-4);
        writetable(Xtbl, [base '_X.csv']);
        writetable(Ytbl, [base '_Y.csv']);
    end
    writetable(table((1:nx)', data.x(:), 'VariableNames', {'index','x'}), ...
               [cfg.filePrefix '_grid_x.csv']);
end

fprintf('\nNext: step3_train_model\n');
