% RUN_SHOCK_TUBE  Solve and plot a single shock-tube problem.
%
% Edit the parameters below (defaults = Sod problem), or set
% interactive = true to be prompted for values like the original script.

clc; clear;

interactive = false;
animate     = false;       % true: animate from t = dt to tEval

% ---- parameters -----------------------------------------------------------
gamma = 1.4;
WL    = [1.0,   0.0, 1.0];     % [rho, u, p] left
WR    = [0.125, 0.0, 0.1];     % [rho, u, p] right
tEval = 0.2;
L     = 1;                     % domain [0, L]
x0    = 0.5;                   % diaphragm position
nx    = 500;
% ---------------------------------------------------------------------------

if interactive
    gamma = input('Gamma: ');
    WL(1) = input('Left density: ');   WL(2) = input('Left velocity: ');
    WL(3) = input('Left pressure: ');
    WR(1) = input('Right density: ');  WR(2) = input('Right velocity: ');
    WR(3) = input('Right pressure: '); tEval  = input('Time of evaluation: ');
end

x = linspace(0, L, nx);

if animate
    times = linspace(tEval / 50, tEval, 50);
else
    times = tEval;
end

labels = {'Density', 'Velocity', 'Pressure'};
for t = times
    [rho, u, p, info] = riemann_exact(gamma, WL, WR, x, t, x0);
    fields = {rho, u, p};
    for k = 1:3
        subplot(1, 3, k);
        plot(x, fields{k}, 'LineWidth', 1.5);
        grid on; xlabel('x'); title(labels{k});
    end
    drawnow;
    if animate, pause(0.02); end
end

fprintf('p* = %.5f,  u* = %.5f\n', info.pStar, info.uStar);
