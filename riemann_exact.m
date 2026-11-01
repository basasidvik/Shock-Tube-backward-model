function [rho, u, p, info] = riemann_exact(gamma, WL, WR, x, t, x0)
%RIEMANN_EXACT  Exact solution of the 1-D Euler shock-tube (Riemann) problem.
%
%   [rho, u, p, info] = riemann_exact(gamma, WL, WR, x, t, x0)
%
%   gamma : ratio of specific heats
%   WL,WR : primitive states [rho, u, p] left / right of the diaphragm
%   x     : array of positions
%   t     : evaluation time (t = 0 returns the initial condition)
%   x0    : diaphragm position (default 0)
%
%   info  : struct with pStar, uStar and the speeds of the outermost left and
%           right waves (sLeft, sRight) - useful to check that the waves are
%           still inside the computational domain.
%
%   Left half is solved directly; the right half is solved by mirroring
%   (x -> -x, u -> -u), so there is only one set of wave formulas.

    if nargin < 6, x0 = 0; end
    WL = WL(:).';  WR = WR(:).';
    if any([WL(1) WL(3) WR(1) WR(3)] <= 0)
        error('riemann:badState', 'Density and pressure must be positive.');
    end

    [pStar, uStar] = riemann_star(gamma, WL, WR);

    WRm = [WR(1), -WR(2), WR(3)];               % mirrored right state
    info.pStar  = pStar;
    info.uStar  = uStar;
    info.sLeft  =  left_edge(gamma, WL,  pStar);
    info.sRight = -left_edge(gamma, WRm, pStar);

    rho = zeros(size(x));  u = rho;  p = rho;

    if t <= 0                                   % initial condition
        isL = x < x0;
        rho(isL) = WL(1);  u(isL) = WL(2);  p(isL) = WL(3);
        rho(~isL) = WR(1); u(~isL) = WR(2); p(~isL) = WR(3);
        return
    end

    S    = (x - x0) / t;                        % similarity variable x/t
    isL  = S <= uStar;                          % left of the contact

    [rho(isL), u(isL), p(isL)] = sample_left(gamma, WL, pStar, uStar, S(isL));

    [r, um, pp] = sample_left(gamma, WRm, pStar, -uStar, -S(~isL));
    rho(~isL) = r;  u(~isL) = -um;  p(~isL) = pp;
end

% =========================================================================
function [rho, u, p] = sample_left(gamma, WK, pStar, uStar, S)
% Solution for S = x/t <= uStar, given the left state WK.
    gm1 = gamma - 1;  gp1 = gamma + 1;
    rhoK = WK(1);  uK = WK(2);  pK = WK(3);
    aK = sqrt(gamma * pK / rhoK);

    rho = zeros(size(S));  u = rho;  p = rho;

    if pStar > pK                               % ---- left shock ----
        sShock  = left_edge(gamma, WK, pStar);
        rhoStar = rhoK * (pStar / pK + gm1 / gp1) / (gm1 / gp1 * pStar / pK + 1);

        pre = S < sShock;                       % undisturbed state
        rho(pre)  = rhoK;     u(pre)  = uK;     p(pre)  = pK;
        rho(~pre) = rhoStar;  u(~pre) = uStar;  p(~pre) = pStar;

    else                                        % ---- left rarefaction ----
        aStar   = aK * (pStar / pK)^(gm1 / (2 * gamma));
        sHead   = uK - aK;
        sTail   = uStar - aStar;
        rhoStar = rhoK * (pStar / pK)^(1 / gamma);

        pre  = S <= sHead;
        fan  = S > sHead & S < sTail;
        post = S >= sTail;

        rho(pre) = rhoK;  u(pre) = uK;  p(pre) = pK;

        c = 2 / gp1 + gm1 / (gp1 * aK) * (uK - S(fan));
        rho(fan) = rhoK * c.^(2 / gm1);
        u(fan)   = 2 / gp1 * (aK + gm1 / 2 * uK + S(fan));
        p(fan)   = pK * c.^(2 * gamma / gm1);

        rho(post) = rhoStar;  u(post) = uStar;  p(post) = pStar;
    end
end

% =========================================================================
function s = left_edge(gamma, WK, pStar)
% Speed of the outermost (leftmost) wave of the left family.
    rhoK = WK(1);  uK = WK(2);  pK = WK(3);
    aK = sqrt(gamma * pK / rhoK);
    if pStar > pK                               % shock speed
        s = uK - aK * sqrt((gamma + 1) / (2 * gamma) * pStar / pK ...
                           + (gamma - 1) / (2 * gamma));
    else                                        % rarefaction head
        s = uK - aK;
    end
end
