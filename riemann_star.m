function [pStar, uStar] = riemann_star(gamma, WL, WR, tol, maxIter)
%RIEMANN_STAR  Pressure and velocity in the star region of the Riemann problem.
%
%   [pStar, uStar] = riemann_star(gamma, WL, WR)
%
%   WL, WR : primitive states [rho, u, p] left / right of the diaphragm.
%   tol    : relative convergence tolerance (default 1e-8)
%   maxIter: maximum Newton iterations      (default 100)
%
%   Newton-Raphson on  f_L(p) + f_R(p) + (uR - uL) = 0, started from a
%   two-shock approximation based on the primitive-variable guess.

    if nargin < 4, tol     = 1e-8; end
    if nargin < 5, maxIter = 100;  end

    rhoL = WL(1); uL = WL(2); pL = WL(3);
    rhoR = WR(1); uR = WR(2); pR = WR(3);

    aL = sqrt(gamma * pL / rhoL);
    aR = sqrt(gamma * pR / rhoR);

    % Pressure-positivity (vacuum) condition
    if 2 * (aL + aR) / (gamma - 1) <= (uR - uL)
        error('riemann:vacuum', ...
              'Initial states generate a vacuum; no star-region solution.');
    end

    pFloor = 1e-8 * min(pL, pR);                % keeps p positive in Newton

    % --- initial guess --------------------------------------------------
    pPV = 0.5 * (pL + pR) - 0.125 * (uR - uL) * (rhoL + rhoR) * (aL + aR);
    pPV = max(pFloor, pPV);

    AL = 2 / ((gamma + 1) * rhoL);  BL = (gamma - 1) / (gamma + 1) * pL;
    AR = 2 / ((gamma + 1) * rhoR);  BR = (gamma - 1) / (gamma + 1) * pR;
    gL = sqrt(AL / (pPV + BL));
    gR = sqrt(AR / (pPV + BR));
    p  = max(pFloor, (gL * pL + gR * pR - (uR - uL)) / (gL + gR));

    % --- Newton-Raphson -------------------------------------------------
    change = inf;
    for iter = 1:maxIter                        %#ok<NASGU>
        [fL, dfL] = pressure_function(p, gamma, rhoL, pL, aL);
        [fR, dfR] = pressure_function(p, gamma, rhoR, pR, aR);

        pNew   = max(pFloor, p - (fL + fR + uR - uL) / (dfL + dfR));
        change = 2 * abs(pNew - p) / (pNew + p);
        p      = pNew;

        if change < tol, break; end
    end
    if change >= tol
        warning('riemann:noConvergence', ...
                'Newton iteration did not converge (change = %.2e).', change);
    end
    pStar = p;

    % --- star velocity (evaluated at the converged pressure) ------------
    fL = pressure_function(pStar, gamma, rhoL, pL, aL);
    fR = pressure_function(pStar, gamma, rhoR, pR, aR);
    uStar = 0.5 * (uL + uR) + 0.5 * (fR - fL);
end
