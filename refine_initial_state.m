function [WL, WR, info] = refine_initial_state(xSensor, Xmeas, gamma, tEval, x0, theta0, opts)
%REFINE_INITIAL_STATE  Pin-point the initial states by fitting the exact solver.
%
%   [WL, WR, info] = refine_initial_state(xSensor, Xmeas, gamma, tEval, x0, theta0)
%   [WL, WR, info] = refine_initial_state(..., opts)
%
%   Starting from the neural-network estimate theta0, this adjusts the six
%   initial values until riemann_exact reproduces the measured sensor
%   readings at ALL measurement times (nonlinear least squares). Using two
%   or more times gives the fit more information than a single snapshot.
%
%   xSensor : 1 x ns sensor positions (use the real ones, not grid-snapped)
%   Xmeas   : 1 x 3*ns*nt raw readings, one block per time (same order as the
%             alltimes CSV):  [rho(1..ns), u(1..ns), p(1..ns)] at tEval(1),
%             then the same block at tEval(2), ...
%   gamma   : ratio of specific heats (known for your gas)
%   tEval   : 1 x nt measurement times;  x0 : diaphragm position
%   theta0  : 1 x 6 NN estimate [rhoL uL pL rhoR uR pR]
%
%   opts (defaults in brackets)
%       method  ['auto']   'lsqnonlin' (Optimization Toolbox), 'fminsearch'
%                          (base MATLAB); auto = lsqnonlin if available
%       fixU    [[]]       [uL uR] to hold the velocities fixed, e.g. [0 0] for
%                          gas initially at rest (then only 4 unknowns remain)
%       nStarts [5]        starting points (the NN estimate + perturbed copies);
%                          max 5
%       perturb [0.15]     size of the perturbations (about 15 %)
%       maxIter [400]
%       priorWeight [0.01] after the fit, slide along the set of equally good
%                          solutions to the one closest to the NN estimate
%                          (0 = off). Only matters when the sensors alone do
%                          not determine the answer.
%
%   Outputs
%       WL, WR : refined [rho u p] left and right states
%       info   : cost (final sum of squared residuals), cost0 (at theta0),
%                method, ambiguity, nNear, costPrior.
%                ambiguity ~ 0 means all good starts agree. A large value
%                means several different initial states fit the readings
%                equally well, i.e. the sensors cannot pin them down.
%
%   Residuals are log(rho_model/rho_meas), log(p_model/p_meas) and
%   (u_model - u_meas)/c, where c is a typical sound speed.

    if nargin < 7, opts = struct(); end
    d = struct('method', 'auto', 'fixU', [], 'nStarts', 5, ...
               'perturb', 0.15, 'maxIter', 400, 'priorWeight', 0.01);
    nm = fieldnames(opts);
    for k = 1:numel(nm), d.(nm{k}) = opts.(nm{k}); end
    if d.nStarts > 5, error('refine_initial_state:nStarts', 'nStarts is at most 5.'); end

    xSensor = xSensor(:).';  ns = numel(xSensor);
    tEval   = tEval(:).';    nt = numel(tEval);
    Xmeas   = Xmeas(:).';
    if numel(Xmeas) ~= 3 * ns * nt
        error('refine_initial_state:badInput', ...
              'Xmeas must have 3*ns*nt = %d values, got %d.', 3 * ns * nt, numel(Xmeas));
    end
    M    = reshape(Xmeas, 3 * ns, nt).';          % nt x 3*ns, one row per time
    rhoM = M(:, 1:ns);  uM = M(:, ns+1:2*ns);  pM = M(:, 2*ns+1:3*ns);
    if any(rhoM(:) <= 0) || any(pM(:) <= 0)
        error('refine_initial_state:badInput', ...
              'Density and pressure readings must be positive.');
    end

    P.x = xSensor;  P.t = tEval;  P.x0 = x0;  P.gamma = gamma;
    P.rhoM = rhoM;  P.uM = uM;  P.pM = pM;
    P.uScale = median(sqrt(gamma * pM(:) ./ rhoM(:)));

    % free parameters: q = [log rhoL, uL/c, log pL, log rhoR, uR/c, log pR]
    P.free   = true(1, 6);
    P.qfixed = zeros(1, 6);
    if ~isempty(d.fixU)
        P.free([2 5])   = false;
        P.qfixed([2 5]) = d.fixU(:).' / P.uScale;
    end
    q0full = [log(theta0(1)), theta0(2) / P.uScale, log(theta0(3)), ...
              log(theta0(4)), theta0(5) / P.uScale, log(theta0(6))];
    q0 = q0full(P.free);

    patterns = [ 0  0  0  0  0  0;
                 1 -1  1 -1  1 -1;
                -1  1 -1  1 -1  1;
                 1  1 -1 -1  1  1;
                -1 -1  1  1 -1 -1];

    method = d.method;
    if strcmp(method, 'auto')
        if exist('lsqnonlin', 'file'), method = 'lsqnonlin';
        else,                          method = 'fminsearch'; end
    end

    cost0 = sum(residual(q0, P).^2);

    qAll = zeros(d.nStarts, numel(q0));
    costs = zeros(d.nStarts, 1);
    for s = 1:d.nStarts
        qs = q0 + d.perturb * patterns(s, P.free);
        [qAll(s, :), costs(s), used] = run_optimizer(qs, @(z) residual(z, P), ...
                                                     method, d.maxIter);
        method = used;                      % stay on the fallback if one was needed
    end

    [best, ib] = min(costs);
    qBest = qAll(ib, :);

    % optional: among equally good solutions prefer the one nearest the NN
    costPrior = best;
    if d.priorWeight > 0
        w  = sqrt(d.priorWeight);
        fp = @(z) [residual(z, P), w * (z - q0)];
        [qp, ~] = run_optimizer(qBest, fp, method, d.maxIter);
        cp = sum(residual(qp, P).^2);
        if cp <= max(2 * best, best + 1e-3)       % still fits the data
            qBest = qp;  costPrior = cp;
        end
    end

    q = P.qfixed;  q(P.free) = qBest;
    WL = [exp(q(1)), q(2) * P.uScale, exp(q(3))];
    WR = [exp(q(4)), q(5) * P.uScale, exp(q(6))];

    near = costs <= max(2 * best, best + 1e-6);
    info.cost      = best;
    info.cost0     = cost0;
    info.costPrior = costPrior;
    info.method    = method;
    info.nNear     = sum(near);
    info.ambiguity = max(max(qAll(near, :), [], 1) - min(qAll(near, :), [], 1));
end

% -------------------------------------------------------------------------
function [q, cost, used] = run_optimizer(q, fres, method, maxIter)
    used = method;
    if strcmp(method, 'lsqnonlin')
        try
            o = optimoptions('lsqnonlin', 'Display', 'off', ...
                'MaxIterations', maxIter, 'FunctionTolerance', 1e-14, ...
                'StepTolerance', 1e-12);
            [q, cost] = lsqnonlin(fres, q, [], [], o);
            return
        catch
            used = 'fminsearch';            % toolbox/licence missing -> fall back
        end
    end
    o = optimset('Display', 'off', 'MaxIter', 10 * maxIter, ...
                 'MaxFunEvals', 20 * maxIter, 'TolX', 1e-10, 'TolFun', 1e-14);
    f = @(z) sum(fres(z).^2);
    for pass = 1:2                          % restart once: Nelder-Mead can stall
        [q, cost] = fminsearch(f, q, o);
    end
end

% -------------------------------------------------------------------------
function res = residual(qf, P)
    ns = numel(P.x);  nt = numel(P.t);
    big = 1e3 * ones(1, 3 * ns * nt);
    q = P.qfixed;  q(P.free) = qf;
    if any(~isfinite(q)) || any(abs(q([1 3 4 6])) > 30)
        res = big;  return
    end
    WL = [exp(q(1)), q(2) * P.uScale, exp(q(3))];
    WR = [exp(q(4)), q(5) * P.uScale, exp(q(6))];
    res = zeros(1, 3 * ns * nt);
    for it = 1:nt
        try
            [r, u, p] = riemann_exact(P.gamma, WL, WR, P.x, P.t(it), P.x0);
        catch err
            if strcmp(err.identifier, 'riemann:vacuum'), res = big;  return; end
            rethrow(err);
        end
        res((it-1)*3*ns + (1:3*ns)) = [log(r) - log(P.rhoM(it, :)), ...
                                       (u - P.uM(it, :)) / P.uScale, ...
                                       log(p) - log(P.pM(it, :))];
    end
    if any(~isfinite(res)), res = big; end
end
