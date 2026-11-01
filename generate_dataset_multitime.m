function files = generate_dataset_multitime(nSamples, times, opts)
%GENERATE_DATASET_MULTITIME  One full dataset file (all grid points) per time.
%
%   files = generate_dataset_multitime(nSamples, times, opts)
%
%   The SAME random initial states are used for every time, so row k of Y is
%   identical in all files and row k of X is the same case at different times.
%   Cases are filtered once with the largest time (the strictest requirement):
%   a case is kept only if its outermost waves are still inside the tube at
%   max(times), so it is valid at every earlier time too.
%
%   Each file holds a struct "data":
%       data.X   nSamples x 3*nx   [rho(1..nx), u(1..nx), p(1..nx)]  (ALL grid points)
%       data.Y   nSamples x 7      [rhoL uL pL rhoR uR pR gamma]
%       data.x   grid positions;  data.opts.tEval = time of this file
%       data.idxTrain / idxVal / idxTest, data.inputChannels, data.targetNames
%
%   opts fields (defaults in brackets): seed [0], nx [256], L [1], x0 [0.5],
%   gammaRange [1.4 1.4], rhoRange [1e-3 2], pRange [1e-3 2], uRange [1e-3 1e-3],
%   margin [0.05], split [0.8 0.1 0.1], filePrefix ['shock_tube_dataset'],
%   precision ['single']  (class of data.X: 'single' or 'double'),
%   uMode ['uniform'] ('uniform' or 'signedlog'), minVal [1e-3] (|u| >= minVal)

    if nargin < 3, opts = struct(); end
    d = struct('seed', 0, 'nx', 256, 'L', 1, 'x0', 0.5, ...
               'gammaRange', [1.4 1.4], 'rhoRange', [1e-3 2], ...
               'pRange', [1e-3 2], 'uRange', [1e-3 1e-3], 'margin', 0.05, ...
               'split', [0.8 0.1 0.1], 'filePrefix', 'shock_tube_dataset', ...
               'precision', 'single', 'uMode', 'uniform', 'minVal', 1e-3);
    names = fieldnames(opts);
    for k = 1:numel(names), d.(names{k}) = opts.(names{k}); end

    rng(d.seed);
    times = times(:).';
    tMax  = max(times);
    x     = linspace(0, d.L, d.nx);
    xLo   = d.margin * d.L;
    xHi   = d.L - d.margin * d.L;

    % ---- 1. draw valid initial states (filtered at the largest time) ------
    % Only the star state and the outer wave speeds are needed here, so the
    % full solver is not called. The star state is kept for step 3.
    fprintf('Sampling %d valid cases (acceptance is typically 10-20%%) ...\n', nSamples);
    tStart = tic;
    Y = zeros(nSamples, 7);
    pStarAll = zeros(nSamples, 1);  uStarAll = zeros(nSamples, 1);
    count = 0;  attempts = 0;  maxAttempts = 200 * nSamples;
    while count < nSamples
        attempts = attempts + 1;
        if attempts > maxAttempts
            error('generate_dataset_multitime:tooManyRejections', ...
                  'Too many rejected samples. Reduce max(times) or the ranges.');
        end
        gamma = uniform(d.gammaRange);
        WL = [log_uniform(d.rhoRange), sample_u(d.uRange, d.uMode, d.minVal), log_uniform(d.pRange)];
        WR = [log_uniform(d.rhoRange), sample_u(d.uRange, d.uMode, d.minVal), log_uniform(d.pRange)];

        aL = sqrt(gamma * WL(3) / WL(1));  aR = sqrt(gamma * WR(3) / WR(1));
        if 2 * (aL + aR) / (gamma - 1) <= WR(2) - WL(2), continue; end   % vacuum

        [pStar, uStar] = riemann_star(gamma, WL, WR);
        sLeft  =  edge_speed(gamma, WL(1),  WL(2), WL(3), pStar);
        sRight = -edge_speed(gamma, WR(1), -WR(2), WR(3), pStar);
        if d.x0 + sLeft  * tMax < xLo, continue; end
        if d.x0 + sRight * tMax > xHi, continue; end

        count = count + 1;
        Y(count, :) = [WL, WR, gamma];
        pStarAll(count) = pStar;  uStarAll(count) = uStar;
        if mod(count, 5000) == 0
            fprintf('  %d / %d cases accepted (%d attempts, %.0f s)\n', ...
                    count, nSamples, attempts, toc(tStart));
        end
    end
    acceptance = nSamples / attempts;
    fprintf('Sampling done in %.0f s (acceptance %.1f%%)\n', toc(tStart), 100 * acceptance);

    % ---- 2. common train/val/test split for all files ---------------------
    idx    = randperm(nSamples);
    nTrain = round(d.split(1) * nSamples);
    nVal   = round(d.split(2) * nSamples);
    idxTrain = idx(1:nTrain);
    idxVal   = idx(nTrain + 1 : nTrain + nVal);
    idxTest  = idx(nTrain + nVal + 1 : end);

    % ---- 3. evaluate every case at every time, one file per time ----------
    files = cell(numel(times), 1);
    chunk = 5000;                                   % cases per vectorised block
    tStart = tic;
    for it = 1:numel(times)
        t = times(it);
        X = zeros(nSamples, 3 * d.nx, d.precision);
        for c0 = 1:chunk:nSamples
            ii = c0:min(c0 + chunk - 1, nSamples);
            [rho, u, p] = riemann_sample_batch(Y(ii,7), Y(ii,1:3), Y(ii,4:6), ...
                                               pStarAll(ii), uStarAll(ii), x, t, d.x0);
            X(ii, :) = [rho, u, p];
        end

        data = struct();
        data.X = X;
        data.Y = Y;
        data.pStar = pStarAll;
        data.uStar = uStarAll;
        data.x = x;
        data.inputChannels = {'rho', 'u', 'p'};
        data.targetNames   = {'rhoL','uL','pL','rhoR','uR','pR','gamma'};
        data.idxTrain = idxTrain;
        data.idxVal   = idxVal;
        data.idxTest  = idxTest;
        data.opts     = d;
        data.opts.tEval = t;
        data.acceptanceRate = acceptance;

        files{it} = sprintf('%s_t%.3f.mat', d.filePrefix, t);
        save(files{it}, 'data', '-v7');
        el = toc(tStart);
        fprintf('t = %.3f -> %s  (X: %d x %d, %s)  [%d/%d, %.0f s elapsed, ~%.0f s left]\n', ...
                t, files{it}, size(X, 1), size(X, 2), d.precision, it, numel(times), ...
                el, el / it * (numel(times) - it));
    end
    fprintf('Done: %d files, %d cases each (acceptance %.1f%%)\n', ...
            numel(times), nSamples, 100 * acceptance);
end

% -------------------------------------------------------------------------
function v = uniform(r)
    v = r(1) + (r(2) - r(1)) * rand;
end

function v = log_uniform(r)
    v = exp(log(r(1)) + (log(r(2)) - log(r(1))) * rand);
end

function v = sample_u(range, mode, minVal)
% Velocity sampler: never exactly zero (|v| >= minVal).
    switch mode
        case 'uniform'
            v = range(1) + (range(2) - range(1)) * rand;
            if abs(v) < minVal
                if v < 0, v = -minVal; else, v = minVal; end
            end
        case 'signedlog'
            hi  = max(abs(range));
            mag = exp(log(minVal) + (log(hi) - log(minVal)) * rand);
            if rand < 0.5, v = -mag; else, v = mag; end
        otherwise
            error('generate_dataset_multitime:badMode', ...
                  'uMode must be ''uniform'' or ''signedlog''.');
    end
end

function s = edge_speed(gamma, rhoK, uK, pK, pStar)
% Speed of the outermost wave of the left family (same as in riemann_exact).
    aK = sqrt(gamma * pK / rhoK);
    if pStar > pK
        s = uK - aK * sqrt((gamma + 1) / (2 * gamma) * pStar / pK + (gamma - 1) / (2 * gamma));
    else
        s = uK - aK;
    end
end
