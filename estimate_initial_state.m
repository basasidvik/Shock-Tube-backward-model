function R = estimate_initial_state(xSensor, Rdata, gamma, times, x0, opts)
%ESTIMATE_INITIAL_STATE  Sensor readings at the measurement times -> initial states.
%
%   R = estimate_initial_state(xSensor, Rdata, gamma, times, x0)
%   R = estimate_initial_state(xSensor, Rdata, gamma, times, x0, opts)
%
%   xSensor : 1 x ns sensor positions (increasing)
%   Rdata   : ns x 3*nt readings. Columns 1:3 = [rho u p] of every sensor at
%             times(1), columns 4:6 at times(2), ...
%   gamma   : ratio of specific heats;  times : 1 x nt;  x0 : diaphragm position
%
%   opts (defaults in brackets)
%       useNN     [true]   start from the neural-network estimate if the model,
%                          sensor layout and times match the training set
%       modelFile ['sensor_model.mat']
%       fixU      [[]]     [uL uR] to hold the velocities fixed (gas at rest)
%       nStarts   [5],  method ['auto'],  ambTol [0.05]
%
%   R fields: theta0 (start), startSource, note, theta (refined [rhoL uL pL
%   rhoR uR pR]), WL, WR, cost, ambiguity, verdict, method, Xmeas.
%
%   This is the function behind shock_tube_app; it has no graphics.

    if nargin < 6, opts = struct(); end
    d = struct('useNN', true, 'modelFile', 'sensor_model.mat', 'fixU', [], ...
               'nStarts', 5, 'method', 'auto', 'ambTol', 0.05);
    nm = fieldnames(opts);
    for k = 1:numel(nm), d.(nm{k}) = opts.(nm{k}); end

    xSensor = xSensor(:).';  ns = numel(xSensor);
    times   = times(:).';    nt = numel(times);
    if size(Rdata, 1) ~= ns || size(Rdata, 2) ~= 3 * nt
        error('estimate_initial_state:badInput', ...
              'Rdata must be %d sensors x %d columns (rho,u,p per time), got %d x %d.', ...
              ns, 3 * nt, size(Rdata, 1), size(Rdata, 2));
    end
    if any(~isfinite(Rdata(:))) || any(diff(xSensor) <= 0)
        error('estimate_initial_state:badInput', ...
              'Readings must be finite and sensor positions increasing.');
    end
    if any(times <= 0)
        error('estimate_initial_state:badInput', 'Measurement times must be positive.');
    end
    left  = find(xSensor < x0, 1, 'first');         % outermost left sensor
    right = find(xSensor > x0, 1, 'last');          % outermost right sensor
    if isempty(left) || isempty(right)
        error('estimate_initial_state:badInput', ...
              'Need at least one sensor on each side of the diaphragm x0.');
    end
    rhoCols = 3 * (0:nt-1) + 1;  pCols = 3 * (0:nt-1) + 3;
    if any(any(Rdata(:, rhoCols) <= 0)) || any(any(Rdata(:, pCols) <= 0))
        error('estimate_initial_state:badInput', 'Density and pressure readings must be positive.');
    end

    % ---- stack readings in the order used for training ---------------------
    Xmeas = zeros(1, 3 * ns * nt);
    for it = 1:nt
        B = Rdata(:, 3 * (it - 1) + (1:3));
        Xmeas((it - 1) * 3 * ns + (1:3 * ns)) = [B(:,1).', B(:,2).', B(:,3).'];
    end

    % ---- starting point ----------------------------------------------------
    theta0 = [];  note = '';  source = '';
    if d.useNN
        [theta0, note] = nn_start(d.modelFile, Xmeas, xSensor, times);
        if ~isempty(theta0), source = 'neural network'; end
    end
    if isempty(theta0)                              % outermost sensors, earliest time
        B1 = Rdata(:, 1:3);
        theta0 = [B1(left, :), B1(right, :)];
        source = 'sensor readings';
    end
    if ~isempty(d.fixU), theta0([2 5]) = d.fixU(:).'; end

    % ---- physics refinement ------------------------------------------------
    ropts = struct('method', d.method, 'nStarts', d.nStarts);
    if ~isempty(d.fixU), ropts.fixU = d.fixU; end
    [WL, WR, info] = refine_initial_state(xSensor, Xmeas, gamma, times, x0, theta0, ropts);

    R.theta0 = theta0;  R.startSource = source;  R.note = note;
    R.theta = [WL, WR];  R.WL = WL;  R.WR = WR;
    R.cost = info.cost;  R.ambiguity = info.ambiguity;  R.method = info.method;
    R.Xmeas = Xmeas;
    if info.ambiguity < d.ambTol
        R.verdict = 'well determined by the sensors';
    else
        R.verdict = 'AMBIGUOUS: several initial states fit these readings; treat with caution';
    end
end

% -------------------------------------------------------------------------
function [theta0, note] = nn_start(modelFile, Xmeas, xSensor, times)
% Neural-network estimate, only if it was trained for exactly this set-up.
    theta0 = [];  note = '';
    if exist(modelFile, 'file') ~= 2
        note = 'No trained model found; started from the sensor readings.';  return
    end
    try
        M  = load(modelFile, 'stats');
        fn = M.stats.featNames;
        okTimes = true;
        for it = 1:numel(times)
            okTimes = okTimes && any(~cellfun(@isempty, strfind(fn, sprintf('_t%.3f', times(it)))));
        end
        cfg = shock_tube_config();
        okPos = numel(cfg.sensorX) == numel(xSensor) && all(abs(cfg.sensorX(:) - xSensor(:)) < 1e-6);
        if numel(fn) ~= numel(Xmeas) || ~okTimes || ~okPos
            note = ['The trained model was built for other times or sensor positions; ' ...
                    'started from the sensor readings.'];
            return
        end
        th = predict_initial_state(Xmeas, modelFile);
        if all(isfinite(th)) && all(th([1 3 4 6]) > 0)
            theta0 = th(:).';
        else
            note = 'The network estimate was not valid; started from the sensor readings.';
        end
    catch err
        note = ['Neural network not used (' err.message '); started from the sensor readings.'];
    end
end