function [rho, u, p] = riemann_sample_batch(gamma, WL, WR, pStar, uStar, x, t, x0)
%RIEMANN_SAMPLE_BATCH  Exact Riemann solution for MANY cases at once.
%
%   [rho, u, p] = riemann_sample_batch(gamma, WL, WR, pStar, uStar, x, t, x0)
%
%   gamma        : N x 1 (or scalar)
%   WL, WR       : N x 3 initial states [rho u p], one row per case
%   pStar, uStar : N x 1 star-region pressure / velocity (from riemann_star)
%   x            : 1 x m positions;   t > 0;   x0 : diaphragm position
%   rho, u, p    : N x m solution arrays
%
%   Same formulas as riemann_exact (one case at a time), but vectorised over
%   the cases, so there is no loop over cases. The star state is passed in
%   because it does not depend on t: it is solved once per case, not once per
%   case and time.

    N = size(WL, 1);
    gamma = gamma(:) .* ones(N, 1);
    pStar = pStar(:);  uStar = uStar(:);
    S = (x(:).' - x0) / t;                          % 1 x m similarity variable
    isLeft = S <= uStar;                            % N x m

    [rL, uL, pL] = sample_left(gamma, WL(:,1),  WL(:,2), WL(:,3), pStar,  uStar,  S);
    [rR, uR, pR] = sample_left(gamma, WR(:,1), -WR(:,2), WR(:,3), pStar, -uStar, -S);

    rho = rR;  u = -uR;  p = pR;                    % right side is the mirrored left
    rho(isLeft) = rL(isLeft);
    u(isLeft)   = uL(isLeft);
    p(isLeft)   = pL(isLeft);
end

% -------------------------------------------------------------------------
function [rho, u, p] = sample_left(gamma, rhoK, uK, pK, pStar, uStar, S)
% Solution for S = x/t <= uStar given the left state; all case data N x 1, S 1 x m.
    gm1 = gamma - 1;  gp1 = gamma + 1;
    aK    = sqrt(gamma .* pK ./ rhoK);
    ratio = pStar ./ pK;
    isShock = pStar > pK;                           % N x 1

    % ---- shock branch ------------------------------------------------------
    sShock = uK - aK .* sqrt(gp1 ./ (2*gamma) .* ratio + gm1 ./ (2*gamma));
    rhoStS = rhoK .* (ratio + gm1 ./ gp1) ./ (gm1 ./ gp1 .* ratio + 1);
    pre = S < sShock;                               % N x m
    rhoS = rhoK .* pre + rhoStS .* ~pre;
    uS   = uK   .* pre + uStar  .* ~pre;
    pS   = pK   .* pre + pStar  .* ~pre;

    % ---- rarefaction branch ------------------------------------------------
    aStar  = aK .* ratio .^ (gm1 ./ (2*gamma));
    sHead  = uK - aK;
    sTail  = uStar - aStar;
    rhoStR = rhoK .* ratio .^ (1 ./ gamma);
    preR  = S <= sHead;
    postR = S >= sTail;
    fan   = ~preR & ~postR;
    c = 2 ./ gp1 + gm1 ./ (gp1 .* aK) .* (uK - S);
    c(~fan) = 1;                                    % keep powers real outside the fan
    c = max(c, 1e-12);
    rhoR = rhoK .* preR + (rhoK .* c .^ (2 ./ gm1)) .* fan + rhoStR .* postR;
    uR   = uK   .* preR + (2 ./ gp1 .* (aK + gm1 ./ 2 .* uK + S)) .* fan + uStar .* postR;
    pR   = pK   .* preR + (pK .* c .^ (2 .* gamma ./ gm1)) .* fan + pStar .* postR;

    % ---- pick the branch per case ------------------------------------------
    rho = rhoR;  u = uR;  p = pR;
    rho(isShock, :) = rhoS(isShock, :);
    u(isShock, :)   = uS(isShock, :);
    p(isShock, :)   = pS(isShock, :);
end
