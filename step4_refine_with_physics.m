% STEP4_REFINE_WITH_PHYSICS  Sharpen the neural-network estimate of the initial
% values by fitting the exact Riemann solver to the sensor readings at all
% measurement times (cfg.useTimes).
%
% Needs the sensor CSVs from step2 and sensor_model.mat from step3.
% For every test case:  NN estimate -> refine_initial_state -> compare to truth.
% Settings: shock_tube_config.m

clc; clear; close all;
cfg    = shock_tube_config();
prefix = cfg.sensorPrefix;
times  = cfg.useTimes;
x0     = cfg.x0;

ropts.method  = 'auto';           % 'lsqnonlin' (Optimization Toolbox) or 'fminsearch'
ropts.nStarts = 5;
ropts.maxIter = 400;
if cfg.gasAtRest
    ropts.fixU = [cfg.minVal cfg.minVal];   % velocities known: only 4 unknowns left
else
    ropts.fixU = [];
end

% ---- load test cases ------------------------------------------------------
if numel(times) > 1
    xFile = [prefix '_alltimes_X.csv'];           % times side by side
else
    xFile = sprintf('%s_t%.3f_X.csv', prefix, times);
end
X = dlmread(xFile, ',', 1, 1);                    % skip header row and case_id
Y = dlmread([prefix '_Y.csv'], ',', 1, 1);        % rhoL uL pL rhoR uR pR gamma

fid = fopen([prefix '_split.csv']);
C = textscan(fid, '%d %s', 'Delimiter', ',', 'HeaderLines', 1);  fclose(fid);
idx = find(strcmp(C{2}, 'test'));
idx = idx(1:min(cfg.nRefine, numel(idx)));

fid = fopen([prefix '_positions.csv']);
C = textscan(fid, '%s %s %f %f %f', 'Delimiter', ',', 'HeaderLines', 1);  fclose(fid);
xs = C{5}(:).';       % sensor positions actually used (for real data, use your own)

if size(X, 2) ~= 3 * numel(xs) * numel(times)
    error('X has %d columns but %d sensors x 3 variables x %d times were expected.', ...
          size(X, 2), numel(xs), numel(times));
end

% ---- stage 1: neural-network estimate -------------------------------------
Theta0 = predict_initial_state(X(idx, :));
Ytrue  = Y(idx, 1:6);

% ---- stage 2: physics refinement ------------------------------------------
n = numel(idx);
Theta1 = zeros(n, 6);  cost = zeros(n, 1);  amb = zeros(n, 1);
for k = 1:n
    [WL, WR, info] = refine_initial_state(xs, X(idx(k), :), Y(idx(k), 7), ...
                                          times, x0, Theta0(k, :), ropts);
    Theta1(k, :) = [WL, WR];
    cost(k) = info.cost;   amb(k) = info.ambiguity;
    if mod(k, 25) == 0, fprintf('refined %d / %d\n', k, n); end
end

% ---- compare --------------------------------------------------------------
names = {'rhoL', 'uL', 'pL', 'rhoR', 'uR', 'pR'};
e0 = Theta0 - Ytrue;   e1 = Theta1 - Ytrue;
rel = @(e, j) 100 * median(abs(e(:, j)) ./ max(abs(Ytrue(:, j)), 1e-3));

fprintf('\n%-6s | %-26s | %-26s\n', 'target', 'NN only', 'NN + physics refinement');
fprintf('%-6s | %10s %14s | %10s %14s\n', '', 'RMSE', 'median rel.err', 'RMSE', 'median rel.err');
for j = 1:6
    if cfg.gasAtRest && any(j == [2 5])      % velocities are fixed: relative error is meaningless
        fprintf('%-6s | %10.4f %14s | %10.4f %14s\n', names{j}, ...
                sqrt(mean(e0(:, j).^2)), 'n/a', sqrt(mean(e1(:, j).^2)), 'n/a');
    else
        fprintf('%-6s | %10.4f %13.1f%% | %10.4f %13.1f%%\n', names{j}, ...
                sqrt(mean(e0(:, j).^2)), rel(e0, j), ...
                sqrt(mean(e1(:, j).^2)), rel(e1, j));
    end
end

well = amb < cfg.ambTol;
fprintf('\nWell determined by the sensors (ambiguity < %.2f): %d of %d cases\n', ...
        cfg.ambTol, sum(well), n);
if any(well)
    fprintf('  median rel. error rhoL/pL/rhoR/pR after refinement: %.1f%% %.1f%% %.1f%% %.1f%%\n', ...
            100 * median(abs(e1(well, [1 3 4 6])) ./ abs(Ytrue(well, [1 3 4 6])), 1));
end
if any(~well)
    fprintf('Ambiguous cases: %d (several initial states fit the readings equally well)\n', ...
            sum(~well));
    fprintf('  median rel. error rhoL/pL/rhoR/pR after refinement: %.1f%% %.1f%% %.1f%% %.1f%%\n', ...
            100 * median(abs(e1(~well, [1 3 4 6])) ./ abs(Ytrue(~well, [1 3 4 6])), 1));
end

figure;
for j = 1:6
    subplot(2, 3, j);
    plot(Ytrue(:, j), Theta0(:, j), '.', 'Color', [0.65 0.65 0.65]); hold on
    plot(Ytrue(:, j), Theta1(:, j), '.', 'Color', [0 0.45 0.74]);
    lim = [min(Ytrue(:, j)) max(Ytrue(:, j))];
    if lim(1) < lim(2), plot(lim, lim, 'r-'); end
    grid on; xlabel('true'); ylabel('estimate'); title(names{j});
    if j == 1, legend('NN', 'refined', 'ideal', 'Location', 'northwest'); end
end

% ---- save results ---------------------------------------------------------
out = [double(idx(:)), Theta0, Theta1, cost, amb];
hdr = ['case_id,' strjoin(strcat(names, '_nn'), ',') ',' ...
       strjoin(strcat(names, '_refined'), ',') ',fit_cost,ambiguity'];
fid = fopen('shock_tube_refined.csv', 'w');  fprintf(fid, '%s\n', hdr);  fclose(fid);
dlmwrite('shock_tube_refined.csv', out, '-append', 'precision', '%.8g');
fprintf('\nSaved shock_tube_refined.csv\n');
