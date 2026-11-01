function cfg = shock_tube_config()
%SHOCK_TUBE_CONFIG  Every setting of the shock-tube pipeline, in ONE place.
%
%   Edit the values below; step1..step4 and check_datasets all read them.
%   All quantities are dimensionless (tube length 1, diaphragm at 0.5).

    % ---- tube and grid ------------------------------------------------------
    cfg.L  = 1;                  % tube length (walls at x = 0 and x = L)
    cfg.x0 = 0.5;                % diaphragm position
    cfg.nx = 256;                % points per stored profile (.mat files)

    % ---- cases --------------------------------------------------------------
    cfg.nSamples   = 50000;      % number of random initial conditions
    cfg.seed       = 0;
    cfg.gammaRange = [1.4 1.4];
    cfg.minVal     = 1e-3;       % smallest density / pressure / velocity used
    cfg.rhoRange   = [cfg.minVal 2];     % density, sampled log-uniformly
    cfg.pRange     = [cfg.minVal 2];     % pressure, sampled log-uniformly
    cfg.gasAtRest  = false;      % true : uL = uR = minVal for every case (constant!)
                                 % false: velocities vary (default)
    cfg.uRangeFree = [-1 1];     % velocity range when gasAtRest = false
    cfg.uMode      = 'uniform';  % 'uniform'   : uniform in uRangeFree
                                 % 'signedlog' : |u| log-uniform from minVal to
                                 %               max|uRangeFree|, random sign
                                 %               (many near-rest cases)
    % In both modes |u| is never smaller than minVal (no exact zero).
    if cfg.gasAtRest
        cfg.uRange = [cfg.minVal cfg.minVal];
    else
        cfg.uRange = cfg.uRangeFree;
    end
    cfg.margin = 0.05;           % waves must stay this far (fraction of L) from walls
    cfg.split  = [0.8 0.1 0.1];  % train / validation / test

    % ---- times --------------------------------------------------------------
    cfg.allTimes = 0.04 * (1:10);   % ALL stored times: 0.04, 0.08, ... , 0.40
    cfg.useTimes = [0.2 0.4];       % the times used for sensors / training
    cfg.storePrecision = 'single';  % 'single' halves the .mat size; 'double' = full
    cfg.filePrefix = 'shock_tube_dataset';

    % ---- sensors (3 left + 3 right of the diaphragm) ------------------------
    cfg.sensorX      = cfg.x0 + [-0.25 -0.15 -0.05,  0.05 0.15 0.25];
    cfg.sensorVars   = {'rho', 'u', 'p'};
    cfg.sensorNoise  = 0;           % relative noise on sensor readings (0.01 = 1 %)
    cfg.sensorMode   = 'exact';     % 'exact': evaluate the solution exactly at the
                                    %          sensor x (no grid snapping)
                                    % 'grid' : take the nearest stored grid point
    cfg.sensorPrefix = 'shock_tube_sensors';
    cfg.exportFullCsv = false;      % true: also dump all 256 points to CSV (huge)

    % ---- training -----------------------------------------------------------
    cfg.method = 'trainnet';        % 'trainnet' (Deep Learning) or 'fitrnet' (Statistics)
    cfg.hidden = [128 128 64];
    cfg.maxEpochs   = 200;          % trainnet: maximum epochs
    cfg.valPatience = 10;           % trainnet: stop after this many validation checks
                                    %           (every 50 iterations) without improvement
    cfg.lossTol     = 1e-6;         % fitrnet: LossTolerance (stop when the loss change is below this)

    % ---- physics refinement -------------------------------------------------
    cfg.nRefine = 300;              % test cases to refine
    cfg.ambTol  = 0.05;             % ambiguity above this = not well determined
end
