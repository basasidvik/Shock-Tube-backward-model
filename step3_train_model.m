% STEP3_TRAIN_MODEL  Learn the shock-tube initial states from sensor readings.
%
% Input : sensor CSVs from step2 (both times side by side: 36 inputs)
% Target: [rhoL uL pL rhoR uR pR]   (gamma is constant, so it is not learned;
%         with cfg.gasAtRest the two velocities are constant as well)
%
% Back-ends (cfg.method):
%   'trainnet' : one network with 6 outputs        -> Deep Learning Toolbox
%                (needs R2023b or newer; otherwise falls back to fitrnet)
%   'fitrnet'  : six separate neural nets          -> Statistics and Machine
%                                                     Learning Toolbox (R2021a+)
% Settings: shock_tube_config.m

clc; clear; close all;
rng(0);
cfg    = shock_tube_config();
prefix = cfg.sensorPrefix;
method = cfg.method;
hidden = cfg.hidden;

if numel(cfg.useTimes) > 1
    xFile = [prefix '_alltimes_X.csv'];
else
    xFile = sprintf('%s_t%.3f_X.csv', prefix, cfg.useTimes);
end
fprintf('Training on %s\n', xFile);

if strcmp(method, 'trainnet') && ~exist('trainnet', 'file')
    warning('trainnet needs R2023b or newer; falling back to fitrnet.');
    method = 'fitrnet';
end

% ---- load data (first column of every file is case_id) --------------------
fid = fopen(xFile);  header = strsplit(fgetl(fid), ',');  fclose(fid);
featNames = header(2:end);
X = readmatrix(xFile, 'NumHeaderLines', 1);   X = X(:, 2:end);
Y = readmatrix([prefix '_Y.csv'], 'NumHeaderLines', 1);
Y = Y(:, 2:7);                                % rhoL uL pL rhoR uR pR (drop gamma)
S = readtable([prefix '_split.csv']);
isTrain = strcmp(S.split, 'train');
isVal   = strcmp(S.split, 'val');
isTest  = strcmp(S.split, 'test');
fprintf('Cases: %d train / %d val / %d test, %d input features\n', ...
        sum(isTrain), sum(isVal), sum(isTest), size(X, 2));

% ---- transforms: log of density and pressure, then standardise ------------
isU  = startsWith(featNames, 'u_');           % velocity features stay linear
logT = [true false true  true false true];    % log of rhoL, pL, rhoR, pR

Xt = X;   Xt(:, ~isU) = log(Xt(:, ~isU));
Yt = Y;   Yt(:, logT) = log(Yt(:, logT));

muX = mean(Xt(isTrain, :));   sgX = std(Xt(isTrain, :));   sgX(sgX < 1e-8) = 1;
muY = mean(Yt(isTrain, :));   sgY = std(Yt(isTrain, :));   % statistics from TRAIN only
sgY(sgY < 1e-8) = 1;                          % constant targets (e.g. u) stay finite

Xn = (Xt - muX) ./ sgX;
Yn = (Yt - muY) ./ sgY;

% ---- train ----------------------------------------------------------------
switch method
    case 'trainnet'
        layers = featureInputLayer(size(Xn, 2));
        for h = hidden
            layers = [layers; fullyConnectedLayer(h); reluLayer]; %#ok<AGROW>
        end
        layers = [layers; fullyConnectedLayer(6)];

        options = trainingOptions('adam', ...
            'MaxEpochs', cfg.maxEpochs, 'MiniBatchSize', 256, ...
            'InitialLearnRate', 1e-3, ...
            'LearnRateSchedule', 'piecewise', ...
            'LearnRateDropFactor', 0.5, 'LearnRateDropPeriod', 40, ...
            'Shuffle', 'every-epoch', ...
            'ValidationData', {Xn(isVal, :), Yn(isVal, :)}, ...
            'ValidationFrequency', 50, 'ValidationPatience', cfg.valPatience, ...
            'Plots', 'training-progress', 'Verbose', false);

        net = trainnet(Xn(isTrain, :), Yn(isTrain, :), layers, 'mse', options);

    case 'fitrnet'
        mdl = cell(1, 6);
        for j = 1:6
            mdl{j} = fitrnet(Xn(isTrain, :), Yn(isTrain, j), ...
                'LayerSizes', hidden, 'Activations', 'relu', ...
                'Standardize', false, 'IterationLimit', 1000, ...
                'LossTolerance', cfg.lossTol, ...
                'ValidationData', {Xn(isVal, :), Yn(isVal, j)});
            fprintf('fitrnet: target %d/6 done\n', j);
        end
    otherwise
        error('Unknown method "%s".', method);
end

% ---- save model + preprocessing statistics --------------------------------
stats = struct('muX', muX, 'sgX', sgX, 'muY', muY, 'sgY', sgY, ...
               'isU', isU, 'logT', logT, 'featNames', {featNames}, ...
               'method', method);
if strcmp(method, 'trainnet')
    save('sensor_model.mat', 'net', 'stats');
else
    save('sensor_model.mat', 'mdl', 'stats');
end

% ---- evaluate on the test set (raw readings -> raw initial values) --------
Ypred = predict_initial_state(X(isTest, :), 'sensor_model.mat');
Ytrue = Y(isTest, :);

names = {'rhoL', 'uL', 'pL', 'rhoR', 'uR', 'pR'};
err   = Ypred - Ytrue;
fprintf('\n%-6s %10s %10s %8s %14s\n', 'target', 'RMSE', 'MAE', 'R^2', 'median rel.err');
for j = 1:6
    ssTot = sum((Ytrue(:,j) - mean(Ytrue(:,j))).^2);
    if ssTot < 1e-12                          % constant target: R^2 undefined
        fprintf('%-6s %10.4f %10.4f %8s %14s\n', names{j}, ...
                sqrt(mean(err(:,j).^2)), mean(abs(err(:,j))), 'n/a', 'n/a');
        continue
    end
    r2  = 1 - sum(err(:,j).^2) / ssTot;
    rel = median(abs(err(:,j)) ./ max(abs(Ytrue(:,j)), 1e-3));
    fprintf('%-6s %10.4f %10.4f %8.4f %13.1f%%\n', names{j}, ...
            sqrt(mean(err(:,j).^2)), mean(abs(err(:,j))), r2, 100 * rel);
end

figure;
for j = 1:6
    subplot(2, 3, j);
    plot(Ytrue(:,j), Ypred(:,j), '.', 'MarkerSize', 4); hold on
    lim = [min(Ytrue(:,j)) max(Ytrue(:,j))];
    if lim(1) < lim(2), plot(lim, lim, 'r-', 'LineWidth', 1.2); end
    grid on; xlabel('true'); ylabel('predicted'); title(names{j});
end

fprintf('\nNext: step4_refine_with_physics\n');
