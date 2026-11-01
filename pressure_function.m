function [f, df] = pressure_function(p, gamma, rhoK, pK, aK)
%PRESSURE_FUNCTION  Toro's f_K(p) and its derivative for one side (K = L or R).
%
%   f_K(p) is the velocity change across the left/right nonlinear wave when
%   the star-region pressure is p:
%       p >  pK : shock
%       p <= pK : rarefaction
%
%   Inputs : p      trial star pressure
%            gamma  ratio of specific heats
%            rhoK, pK, aK   density, pressure, sound speed of the side state
%   Outputs: f, df  function value and derivative d f / d p

    if p > pK                                   % shock
        A  = 2 / ((gamma + 1) * rhoK);
        B  = (gamma - 1) / (gamma + 1) * pK;
        q  = sqrt(A / (p + B));
        f  = (p - pK) * q;
        df = q * (1 - (p - pK) / (2 * (B + p)));
    else                                        % rarefaction
        z  = (gamma - 1) / (2 * gamma);
        f  = 2 * aK / (gamma - 1) * ((p / pK)^z - 1);
        df = (p / pK)^(-(gamma + 1) / (2 * gamma)) / (rhoK * aK);
    end
end
