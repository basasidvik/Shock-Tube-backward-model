% INSPECT_MODEL  Show and export the trained network's weights and biases.
%
% Run after step3_train_model (which creates sensor_model.mat). Step 4 does not
% change the model, it only uses it. Prints every layer's size, writes each
% weight matrix / bias vector to a CSV (model_*.csv), and writes the input and
% output normalisation that the weights assume.

clc; clear;
M = load('sensor_model.mat');
s = M.stats;
fprintf('Model type: %s\n', s.method);
fprintf('%d input features: %s ... %s\n', numel(s.featNames), s.featNames{1}, s.featNames{end});

total = 0;
if strcmp(s.method, 'trainnet')
    L = M.net.Learnables;                     % table: Layer, Parameter, Value
    fprintf('\n%-12s %-10s %s\n', 'layer', 'parameter', 'size');
    for i = 1:height(L)
        v = double(extractdata(L.Value{i}));
        total = total + numel(v);
        lname = char(L.Layer(i));  pname = char(L.Parameter(i));
        fprintf('%-12s %-10s %s\n', lname, pname, mat2str(size(v)));
        writematrix(v, sprintf('model_%s_%s.csv', lname, pname));
    end
else                                          % fitrnet: one network per target
    names = {'rhoL', 'uL', 'pL', 'rhoR', 'uR', 'pR'};
    for j = 1:numel(M.mdl)
        W = M.mdl{j}.LayerWeights;  B = M.mdl{j}.LayerBiases;
        fprintf('\nTarget %s: layer sizes %s\n', names{j}, mat2str(M.mdl{j}.LayerSizes));
        for k = 1:numel(W)
            total = total + numel(W{k}) + numel(B{k});
            fprintf('  layer %d: weights %s, biases %s\n', k, ...
                mat2str(size(W{k})), mat2str(size(B{k})));
            writematrix(W{k}, sprintf('model_%s_layer%d_weights.csv', names{j}, k));
            writematrix(B{k}, sprintf('model_%s_layer%d_biases.csv',  names{j}, k));
        end
    end
end
fprintf('\nTotal learnable parameters: %d\n', total);

% ---- normalisation used before / after the network ------------------------
writematrix([s.muX; s.sgX], 'model_input_mean_std.csv');     % row 1 mean, row 2 std
writematrix([s.muY; s.sgY; double(s.logT)], 'model_output_mean_std_logflag.csv');
fprintf('Wrote model_*.csv (weights, biases, normalisation) to %s\n', pwd);