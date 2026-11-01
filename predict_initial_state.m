function Y = predict_initial_state(Xraw, modelFile)
%PREDICT_INITIAL_STATE  Predict shock-tube initial states from sensor readings.
%
%   Y = predict_initial_state(Xraw)
%   Y = predict_initial_state(Xraw, modelFile)
%
%   Xraw : N x F matrix of RAW sensor readings, columns in the same order as
%          the training CSV (rho_S1..rho_S6, u_S1..u_S6, p_S1..p_S6, repeated
%          for each time if the model was trained on the alltimes file).
%          Density and pressure readings must be > 0.
%   Y    : N x 6 matrix  [rhoL uL pL rhoR uR pR]  (same units as the data
%          the model was trained on)
%
%   modelFile defaults to 'sensor_model.mat' (written by train_sensor_model).

    if nargin < 2, modelFile = 'sensor_model.mat'; end
    M = load(modelFile);
    s = M.stats;

    if size(Xraw, 2) ~= numel(s.isU)
        error('predict_initial_state:badInput', ...
              'Expected %d columns, got %d.', numel(s.isU), size(Xraw, 2));
    end
    if any(any(Xraw(:, ~s.isU) <= 0))
        error('predict_initial_state:badInput', ...
              'Density and pressure readings must be positive.');
    end

    X = Xraw;
    X(:, ~s.isU) = log(X(:, ~s.isU));
    Xn = (X - s.muX) ./ s.sgX;

    switch s.method
        case 'trainnet'
            Yn = predict(M.net, Xn);
            if isa(Yn, 'dlarray'), Yn = extractdata(Yn); end
            Yn = double(Yn);
            if size(Yn, 1) ~= size(Xn, 1), Yn = Yn.'; end
        case 'fitrnet'
            Yn = zeros(size(Xn, 1), numel(M.mdl));
            for j = 1:numel(M.mdl)
                Yn(:, j) = predict(M.mdl{j}, Xn);
            end
        otherwise
            error('Unknown method "%s" in model file.', s.method);
    end

    Y = Yn .* s.sgY + s.muY;
    Y(:, s.logT) = exp(Y(:, s.logT));
end
