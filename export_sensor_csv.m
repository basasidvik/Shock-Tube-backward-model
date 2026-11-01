function export_sensor_csv(files, sensorX, opts)
%EXPORT_SENSOR_CSV  Labelled CSVs of the readings of a few "sensors".
%
%   export_sensor_csv(files, sensorX)
%   export_sensor_csv(files, sensorX, opts)
%
%   files   : cell array of the .mat dataset files to use (e.g. the 2 times
%             used for training). The .mat files are NOT modified.
%   sensorX : sensor positions (e.g. 3 left + 3 right of the diaphragm)
%
%   opts (defaults in brackets)
%       mode     ['exact']  'exact': evaluate the exact solution AT the sensor
%                           positions (no grid snapping, no interpolation);
%                           'grid': nearest stored grid point
%       vars     [{'rho','u','p'}]
%       noiseStd [0]        relative Gaussian sensor noise (0.01 = 1 %)
%       seed     [1]
%       prefix   ['shock_tube_sensors']
%
%   Files written (labelled, with a case_id column):
%       <prefix>_t0.200_X.csv ...   readings, one file per time
%       <prefix>_alltimes_X.csv     all times side by side (if > 1 file)
%       <prefix>_Y.csv              targets: rhoL uL pL rhoR uR pR gamma
%       <prefix>_split.csv          train / val / test label per case
%       <prefix>_positions.csv      sensor name, side, grid index, x requested/used
%   Columns are <var>_S<k> (sensors numbered left to right), with a time
%   suffix in the alltimes file, e.g. p_S3_t0.400.

    if nargin < 3, opts = struct(); end
    d = struct('mode', 'exact', 'vars', {{'rho', 'u', 'p'}}, 'noiseStd', 0, ...
               'seed', 1, 'prefix', 'shock_tube_sensors');
    nm = fieldnames(opts);
    for k = 1:numel(nm), d.(nm{k}) = opts.(nm{k}); end
    rng(d.seed);

    chOffset = struct('rho', 0, 'u', 1, 'p', 2);
    for iv = 1:numel(d.vars)
        if ~isfield(chOffset, d.vars{iv})
            error('export_sensor_csv:badVar', 'Unknown variable "%s".', d.vars{iv});
        end
    end
    if ~any(strcmp(d.mode, {'exact', 'grid'}))
        error('export_sensor_csv:badMode', 'mode must be ''exact'' or ''grid''.');
    end

    sensorX = sort(sensorX(:).');
    ns      = numel(sensorX);
    first   = load(files{1});
    x       = first.data.x;
    x0      = first.data.opts.x0;
    nx      = numel(x);
    Y       = double(first.data.Y);
    nCases  = size(Y, 1);

    [~, gi] = min(abs(x(:) - sensorX), [], 1);      % nearest grid points
    if strcmp(d.mode, 'grid')
        if numel(unique(gi)) < ns
            error('export_sensor_csv:duplicateSensors', ...
                  'Two sensors fall on the same grid point.');
        end
        xUse = x(gi);
    else
        xUse = sensorX;
    end
    isLeft = xUse < x0;

    allX = [];  allNames = {};
    for it = 1:numel(files)
        S = load(files{it});  data = S.data;
        if ~isequal(double(data.Y), Y)
            error('export_sensor_csv:mismatch', 'Files hold different cases (Y differs).');
        end
        t = data.opts.tEval;

        % ---- readings: rho, u, p at the sensors -> nCases x ns each ----------
        R.rho = zeros(nCases, ns);  R.u = R.rho;  R.p = R.rho;
        if strcmp(d.mode, 'exact')
            if isfield(data, 'pStar')            % fast: vectorised over all cases
                [R.rho, R.u, R.p] = riemann_sample_batch(Y(:,7), Y(:,1:3), Y(:,4:6), ...
                    double(data.pStar), double(data.uStar), sensorX, t, x0);
            else                                 % older files: one case at a time
                for k = 1:nCases
                    [R.rho(k,:), R.u(k,:), R.p(k,:)] = riemann_exact( ...
                        Y(k,7), Y(k,1:3), Y(k,4:6), sensorX, t, x0);
                end
            end
        else
            R.rho = double(data.X(:, gi));
            R.u   = double(data.X(:, nx + gi));
            R.p   = double(data.X(:, 2*nx + gi));
        end

        % how many cases still show the undisturbed initial state at each sensor
        und = zeros(1, ns);
        for is = 1:ns
            vals = [R.rho(:, is), R.u(:, is), R.p(:, is)];
            if isLeft(is), ref = Y(:, 1:3); else, ref = Y(:, 4:6); end
            und(is) = mean(all(abs(vals - ref) <= 1e-6 * abs(ref) + 1e-12, 2));
        end
        fprintf('t = %.3f: undisturbed fraction per sensor: %s\n', t, sprintf('%.2f ', und));

        % ---- assemble in the documented column order -------------------------
        Xs    = zeros(nCases, numel(d.vars) * ns);
        names = cell(1, numel(d.vars) * ns);
        c = 0;
        for iv = 1:numel(d.vars)
            for is = 1:ns
                c = c + 1;
                Xs(:, c) = R.(d.vars{iv})(:, is);
                names{c} = sprintf('%s_S%d', d.vars{iv}, is);
            end
        end
        if d.noiseStd > 0
            Xs = Xs .* (1 + d.noiseStd * randn(size(Xs)));
        end

        write_numeric_csv(sprintf('%s_t%.3f_X.csv', d.prefix, t), names, Xs);
        allX     = [allX, Xs];                                          %#ok<AGROW>
        allNames = [allNames, cellfun(@(s) sprintf('%s_t%.3f', s, t), ...
                                      names, 'UniformOutput', false)];  %#ok<AGROW>
    end

    if numel(files) > 1
        write_numeric_csv([d.prefix '_alltimes_X.csv'], allNames, allX);
    end
    write_numeric_csv([d.prefix '_Y.csv'], first.data.targetNames, Y);

    split = cell(nCases, 1);
    split(first.data.idxTrain) = {'train'};
    split(first.data.idxVal)   = {'val'};
    split(first.data.idxTest)  = {'test'};
    fid = fopen([d.prefix '_split.csv'], 'w');
    fprintf(fid, 'case_id,split\n');
    for k = 1:nCases, fprintf(fid, '%d,%s\n', k, split{k}); end
    fclose(fid);

    fid = fopen([d.prefix '_positions.csv'], 'w');
    fprintf(fid, 'sensor,side,grid_index,x_requested,x_actual\n');
    for is = 1:ns
        side = 'right';  if isLeft(is), side = 'left'; end
        fprintf(fid, 'S%d,%s,%d,%.6g,%.6g\n', is, side, gi(is), sensorX(is), xUse(is));
    end
    fclose(fid);

    fprintf('Sensor CSVs written (%s mode, %d sensors x %d variables, %d cases, prefix "%s")\n', ...
            d.mode, ns, numel(d.vars), nCases, d.prefix);
end

% -------------------------------------------------------------------------
function write_numeric_csv(file, header, M)
    fid = fopen(file, 'w');
    fprintf(fid, '%s\n', strjoin([{'case_id'}, header], ','));
    fclose(fid);
    dlmwrite(file, [(1:size(M, 1))', M], '-append', 'precision', '%.8g');
end
