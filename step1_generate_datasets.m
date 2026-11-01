% STEP1_GENERATE_DATASETS  Generate and store the FULL datasets (all grid points).
%
% Writes one .mat file per time in cfg.allTimes (10 files: t = 0.04 ... 0.40),
% each with cfg.nSamples cases and cfg.nx grid points x 3 variables per case.
% This is the long step (several minutes for 50,000 cases); run it once.
% Settings: shock_tube_config.m

clc; clear;
cfg = shock_tube_config();

opts = struct('seed', cfg.seed, 'nx', cfg.nx, 'L', cfg.L, 'x0', cfg.x0, ...
              'gammaRange', cfg.gammaRange, 'rhoRange', cfg.rhoRange, ...
              'pRange', cfg.pRange, 'uRange', cfg.uRange, ...
              'margin', cfg.margin, 'split', cfg.split, ...
              'filePrefix', cfg.filePrefix, 'precision', cfg.storePrecision, ...
              'uMode', cfg.uMode, 'minVal', cfg.minVal);

fprintf('Generating %d cases at %d times, %d grid points each ...\n', ...
        cfg.nSamples, numel(cfg.allTimes), cfg.nx);
files = generate_dataset_multitime(cfg.nSamples, cfg.allTimes, opts);

fprintf('\nStored %d full .mat files. Next: step2_export_sensor_csv\n', numel(files));
