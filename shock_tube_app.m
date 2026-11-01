function shock_tube_app()
%SHOCK_TUBE_APP  App: sensor readings at the measurement times -> initial states.
%
%   Run:  shock_tube_app
%
%   Type (or paste) the readings of the sensors at each measurement time into
%   the table, press "Estimate initial states". The app
%     1. picks a starting guess (neural network if its model matches this
%        set-up, otherwise the readings of the outermost sensors),
%     2. refines it by fitting the exact Riemann solution to all readings,
%     3. shows the six initial values and re-plots the flow at each time with
%        your readings overlaid.
%
%   Needs shock_tube_config.m (sensor layout and times) and, optionally,
%   sensor_model.mat from step3. Values are dimensionless unless you use your
%   own consistent units. The exact model is valid only until the first wave
%   reaches a wall; the app warns when the result violates that.

    cfg    = shock_tube_config();
    ns     = numel(cfg.sensorX);
    nt     = numel(cfg.useTimes);
    names  = {'rhoL', 'uL', 'pL', 'rhoR', 'uR', 'pR'};
    trueVals = [];                                    % known truth (example / dataset case)
    haveNN = exist('sensor_model.mat', 'file') == 2;

    % ---------------------------- layout -------------------------------------
    fig  = uifigure('Name', 'Shock tube - initial-state estimator', ...
                    'Position', [40 40 1200 740]);
    root = uigridlayout(fig, [1 2]);
    root.ColumnWidth = {500, '1x'};

    left = uigridlayout(root, [5 1]);
    left.RowHeight = {'fit', '1x', 'fit', 'fit', 'fit'};

    % --- settings
    nc = max(4, nt + 1);
    sg = uigridlayout(left, [5 nc]);
    sg.RowHeight = repmat({'fit'}, 1, 5);
    uilabel(sg, 'Text', 'gamma');
    gammaField = uieditfield(sg, 'numeric', 'Value', 1.4, 'Limits', [1.0001 3]);
    uilabel(sg, 'Text', 'x0 (diaphragm)');
    x0Field = uieditfield(sg, 'numeric', 'Value', cfg.x0);
    lbl = uilabel(sg, 'Text', 'Tube length L');
    LField = uieditfield(sg, 'numeric', 'Value', cfg.L, 'Limits', [1e-9 Inf]);
    lbl2 = uilabel(sg, 'Text', 'Times');
    lbl2.Layout.Row = 3;  lbl2.Layout.Column = 1;
    timeFields = cell(1, nt);
    for it = 1:nt
        timeFields{it} = uieditfield(sg, 'numeric', 'Value', cfg.useTimes(it), ...
                                     'Limits', [1e-12 Inf], 'ValueChangedFcn', @(~,~) updateColNames());
        timeFields{it}.Layout.Row = 3;  timeFields{it}.Layout.Column = 1 + it;
    end
    nnBox = uicheckbox(sg, 'Text', 'Use the neural network as starting point (if it matches)', ...
                       'Value', haveNN, 'Enable', onoff(haveNN));
    nnBox.Layout.Row = 4;  nnBox.Layout.Column = [1 nc];
    restBox = uicheckbox(sg, 'Text', 'Gas initially at rest: fix both velocities to');
    restBox.Layout.Row = 5;  restBox.Layout.Column = [1 nc - 1];
    restField = uieditfield(sg, 'numeric', 'Value', cfg.minVal);
    restField.Layout.Row = 5;  restField.Layout.Column = nc;

    % --- readings table
    rowNames = arrayfun(@(k) sprintf('S%d', k), 1:ns, 'UniformOutput', false);
    tbl = uitable(left, 'Data', [cfg.sensorX(:), zeros(ns, 3 * nt)], ...
                  'RowName', rowNames, 'ColumnEditable', true(1, 1 + 3 * nt));

    % --- buttons
    bg = uigridlayout(left, [1 3]);
    uibutton(bg, 'Text', 'Estimate initial states', 'FontWeight', 'bold', ...
             'ButtonPushedFcn', @onEstimate);
    uibutton(bg, 'Text', 'Load example', 'ButtonPushedFcn', @onExample);
    uibutton(bg, 'Text', 'Clear readings', 'ButtonPushedFcn', @onClear);

    cg = uigridlayout(left, [1 3]);
    cg.ColumnWidth = {'fit', 90, '1x'};
    uilabel(cg, 'Text', 'Dataset case #');
    caseField = uieditfield(cg, 'numeric', 'Value', 1, 'Limits', [1 Inf], ...
                            'RoundFractionalValues', true);
    uibutton(cg, 'Text', 'Load case from sensor CSVs', 'ButtonPushedFcn', @onLoadCase);

    statusLbl = uilabel(left, 'Text', 'Enter the readings or load an example.', 'WordWrap', 'on');

    % --- results and plots
    right = uigridlayout(root, [5 1]);
    right.RowHeight = {150, 'fit', '1x', '1x', '1x'};
    resTbl = uitable(right, 'ColumnName', {'Initial value', 'Start', 'Refined', 'True'}, ...
                     'RowName', {}, 'Data', cell(6, 4));
    summary = uilabel(right, 'Text', '', 'WordWrap', 'on');
    ax = gobjects(1, 3);
    for k = 1:3, ax(k) = uiaxes(right); end
    updateColNames();

    % =========================== callbacks ====================================
    function onEstimate(~, ~)
        try
            setStatus('Working ...');  drawnow;
            [x, Rdata, gamma, tms, x0v, Lt] = readInputs();
            opts = struct('useNN', nnBox.Value);
            if restBox.Value, opts.fixU = [restField.Value restField.Value]; end
            res = estimate_initial_state(x, Rdata, gamma, tms, x0v, opts);
            showResult(res, gamma, tms, x, Rdata, x0v, Lt);
            setStatus('Done.');
        catch err
            setStatus(['Error: ' err.message]);
            uialert(fig, err.message, 'Could not estimate the initial states');
        end
    end

    function onExample(~, ~)
        gamma = 1.4;  WL = [1.0 0.001 0.5];  WR = [0.5 0.001 0.2];   % gentle example case
        tms = cfg.useTimes;  x = cfg.sensorX;
        Rd = zeros(ns, 3 * nt);
        for it = 1:nt
            [r, u, p] = riemann_exact(gamma, WL, WR, x, tms(it), cfg.x0);
            Rd(:, 3 * (it - 1) + (1:3)) = [r(:), u(:), p(:)];
        end
        fillInputs(x, Rd, gamma, tms);
        trueVals = [WL WR];
        setStatus('Example loaded (exact readings from known initial states). Press Estimate.');
    end

    function onClear(~, ~)
        tbl.Data(:, 2:end) = 0;
        trueVals = [];
        resTbl.Data = cell(6, 4);  summary.Text = '';
        for k = 1:3, cla(ax(k)); end
        setStatus('Cleared. Enter the readings.');
    end

    function onLoadCase(~, ~)
        try
            k = round(caseField.Value);
            pre = cfg.sensorPrefix;
            if nt > 1, xFile = [pre '_alltimes_X.csv'];
            else,      xFile = sprintf('%s_t%.3f_X.csv', pre, cfg.useTimes); end
            row = dlmread(xFile, ',', [k 1 k 3 * ns * nt]);
            y   = dlmread([pre '_Y.csv'], ',', [k 1 k 7]);
            Rd  = zeros(ns, 3 * nt);
            for it = 1:nt
                b = row((it - 1) * 3 * ns + (1:3 * ns));
                Rd(:, 3 * (it - 1) + (1:3)) = [b(1:ns).', b(ns+1:2*ns).', b(2*ns+1:3*ns).'];
            end
            fillInputs(cfg.sensorX, Rd, y(7), cfg.useTimes);
            trueVals = y(1:6);
            setStatus(sprintf('Case %d loaded; the true initial values are shown after Estimate.', k));
        catch err
            setStatus(['Error: ' err.message]);
            uialert(fig, ['Could not read the sensor CSVs (run step2 first): ' err.message], 'Load case');
        end
    end

    % =========================== helpers ======================================
    function [x, Rdata, gamma, tms, x0v, Lt] = readInputs()
        D     = tbl.Data;
        x     = D(:, 1).';
        Rdata = D(:, 2:end);
        gamma = gammaField.Value;  x0v = x0Field.Value;  Lt = LField.Value;
        tms   = cellfun(@(f) f.Value, timeFields);
    end

    function fillInputs(x, Rd, gamma, tms)
        tbl.Data = [x(:), Rd];
        gammaField.Value = gamma;  x0Field.Value = cfg.x0;  LField.Value = cfg.L;
        for it = 1:nt, timeFields{it}.Value = tms(it); end
        updateColNames();
        resTbl.Data = cell(6, 4);  summary.Text = '';
        for k = 1:3, cla(ax(k)); end
    end

    function updateColNames()
        cn = {'x'};
        for it = 1:nt
            t = timeFields{it}.Value;
            cn = [cn, {sprintf('rho (t=%g)', t), sprintf('u (t=%g)', t), sprintf('p (t=%g)', t)}]; %#ok<AGROW>
        end
        tbl.ColumnName = cn;
    end

    function setStatus(msg)
        statusLbl.Text = msg;
    end

    function showResult(res, gamma, tms, x, Rdata, x0v, Lt)
        d = cell(6, 4);
        for i = 1:6
            d{i, 1} = names{i};
            d{i, 2} = sprintf('%.5g', res.theta0(i));
            d{i, 3} = sprintf('%.5g', res.theta(i));
            if isempty(trueVals), d{i, 4} = '-'; else, d{i, 4} = sprintf('%.5g', trueVals(i)); end
        end
        resTbl.Data = d;

        txt = sprintf('Start: %s | fit cost %.2e | ambiguity %.3f | %s (%s)', ...
                      res.startSource, res.cost, res.ambiguity, res.verdict, res.method);
        if ~isempty(res.note), txt = [txt newline 'Note: ' res.note]; end
        tmax = max(tms);
        [~, ~, ~, inf0] = riemann_exact(gamma, res.WL, res.WR, x0v, tmax, x0v);
        if x0v + inf0.sLeft * tmax < 0 || x0v + inf0.sRight * tmax > Lt
            txt = [txt newline sprintf(['Warning: at t = %g the refined waves would already have ' ...
                   'reached a wall, so the exact model is not valid there.'], tmax)];
        end
        summary.Text = txt;

        % flow at each time from the refined initial states, readings overlaid
        xs = linspace(0, Lt, 400);
        [r0, u0, p0] = riemann_exact(gamma, res.WL, res.WR, xs, 0, x0v);
        prof0 = {r0, u0, p0};
        prof  = cell(1, nt);
        for it = 1:nt
            [r, u, p] = riemann_exact(gamma, res.WL, res.WR, xs, tms(it), x0v);
            prof{it} = {r, u, p};
        end
        cols = lines(nt);  ylab = {'density', 'velocity', 'pressure'};
        for k = 1:3
            a = ax(k);  cla(a);  hold(a, 'on');
            h0 = plot(a, xs, prof0{k}, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2);
            hs = gobjects(1, nt);
            for it = 1:nt
                hs(it) = plot(a, xs, prof{it}{k}, '-', 'Color', cols(it, :), 'LineWidth', 1.4);
                plot(a, x, Rdata(:, 3 * (it - 1) + k), 'o', 'Color', cols(it, :), ...
                     'MarkerFaceColor', cols(it, :));
            end
            ylabel(a, ylab{k});  xlabel(a, 'x');  grid(a, 'on');  xlim(a, [0 Lt]);
            if k == 1
                legend(a, [h0, hs], [{'t = 0 (estimated)'}, ...
                       arrayfun(@(t) sprintf('t = %g', t), tms, 'UniformOutput', false)], ...
                       'Location', 'best');
            end
            hold(a, 'off');
        end
    end
end

function s = onoff(tf)
    if tf, s = 'on'; else, s = 'off'; end
end