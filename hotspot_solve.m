function sol = hotspot_solve(p, opts)
%HOTSPOT_SOLVE  Tau-driven moving hotspot of Abeta secondary nucleation, 1-D.
%
%   sol = HOTSPOT_SOLVE(p, opts)
%
% Integrates, on 0 <= x <= 1 with no-flux boundaries,
%
%   dT/dt = D_T T_xx + r_T T (1 - T)                                      (11)
%   dm/dt = D_m m_xx + lam (1 - m) - 2 k_E m P - nu (R_tau + R_2)         (12)
%   dP/dt = R_tau + R_2 - mu_P P                                          (13)
%   dM/dt = 2 k_E m P + nu (R_tau + R_2) - mu_M M                         (14)
%   R_tau = k_tau T m^2,   R_2 = k_2 A(M) m^2,   A = M/(1+(M/M_c)^(1/3))  (15)
%
% With opts.qss = true (default) Eq. (12) is replaced by its quasi-steady form
%
%   0 = D_m m_xx + lam (1 - m) - 2 k_E m P - nu (R_tau + R_2)             (16)
%
% which is solved for m by damped Newton iteration at every evaluation of the
% right-hand side. This is justified by the separation of time scales set out
% in Section 2.6: monomer relaxes in ~1 h while the tau field evolves over
% months, a factor of ~1.4e3 in mouse and ~7e4 in human. Slaving the monomer
% removes the stiffest scale and leaves the time unit fixed by the tau field,
% which is what makes the dimensional calibration possible.
%
% opts fields (all optional)
%   N        grid points                              (default 161)
%   tspan    output times                             (default 0:0.1:30)
%   qss      quasi-steady monomer                     (default true)
%   RelTol   ode15s relative tolerance                (default 2e-6)
%   AbsTol   ode15s absolute tolerance                (default 1e-9)
%   newtTol  Newton residual tolerance for Eq. (16)   (default 1e-11)
%   verbose  print progress                           (default false)
%
% Returns sol with fields x, t, T, m, P, M (each n_t by N) and sol.stats.
%
% The monomer operator in Eq. (16) is negative definite for m > 0, lam > 0 and
% non-negative nucleation coefficients, so the positive solution is unique and
% Newton converges from any positive start.

if nargin < 1 || isempty(p), p = hotspot_params('reference'); end
if nargin < 2, opts = struct(); end
d = @(f, v) ternary(isfield(opts, f), @() opts.(f), @() v);
N       = d('N', 161);
tspan   = d('tspan', 0:0.1:30);
qss     = d('qss', true);
RelTol  = d('RelTol', 2e-6);
AbsTol  = d('AbsTol', 1e-9);
newtTol = d('newtTol', 1e-11);
verbose = d('verbose', false);

% ---- grid and Neumann Laplacian -------------------------------------------
x = linspace(0, 1, N).';
h = x(2) - x(1);
e = ones(N, 1);
L = spdiags([e, -2*e, e], -1:1, N, N);
L(1, 2)     = 2;      % ghost node m_0 = m_2  (zero flux)
L(N, N-1)   = 2;
L = L / h^2;

% Optional diffusing oligomer intermediate. Peter et al. (2026) showed that
% spatial transmission of Abeta aggregation is carried by oligomers diffusing
% ahead of the front, not by fibril growth. With oligo = true the nucleation
% flux feeds a mobile species S which converts to fibrils at rate k_conv:
%     dS/dt = D_S S_xx + (R_tau + R_2) - k_conv S
%     dP/dt = k_conv S - mu_P P
%     dM/dt = 2 k_E m P + nu k_conv S - mu_M M
% recovering the base model as D_S -> 0 with k_conv -> Inf. The oligomer
% transport length is sqrt(D_S/k_conv).
oligo = isfield(p, 'D_S') && isfield(p, 'k_conv') && p.D_S > 0;

% Optional history-dependent catalytic-site density. Hu et al. (2026) found
% fibrils annealed at low supersaturation to carry roughly 8-fold fewer
% secondary-nucleation defects than fibrils grown rapidly, with no change in
% morphology. With sig_lo < 1 the site density laid down at each instant is
%     sigma(m) = sig_lo + (1 - sig_lo) m
% normalised to 1 at the background monomer level, and a separate site-weighted
% mass Q accumulates with that weight while M accumulates without it. Exposed
% catalytic capacity is then Q/(1 + (M/M_c)^(1/3)) instead of A(M), so that
% fibrils born in the depleted wake are intrinsically less catalytic.
hist = isfield(p, 'sig_lo') && p.sig_lo < 1;

% Optional reversible deposition. Cirrito et al. (2003) found the apparent
% half-life of interstitial Abeta to be roughly DOUBLED in plaque-bearing mice,
% which an irreversible sink cannot produce: deposits behave as reservoirs that
% return monomer. With k_off > 0 the elongation term becomes a net flux
%     net = 2 k_E m P - k_off M
% so that monomer is replenished from fibril mass as it dissolves, and the
% local equilibrium monomer concentration is k_off M / (2 k_E P).

% ---- initial state ---------------------------------------------------------
T0 = p.T0 * exp(-(x / p.w).^2);
m0 = ones(N, 1);
P0 = zeros(N, 1);
M0 = zeros(N, 1);

mWarm = m0;            % Newton warm start, persists across RHS calls
nNewton = 0; nRHS = 0; maxIter = 0;

    function A = accessible(M)
        % A(M) = M / (1 + (M/M_c)^(1/3)), guarded at M = 0
        Mp = max(M, 0);
        A  = Mp ./ (1 + (Mp / p.M_c).^(1/3));
    end

    function m = solve_monomer(T, P, M, mGuess)
        % Damped Newton on F(m) = D_m L m + lam(1-m) - 2 k_E m P - nu*b*m^2
        b = p.nu * (p.k_tau * T + p.k_2 * accessible(M));
        a = p.lam + 2 * p.k_E * P;
        if isfield(p, 'k_off') && p.k_off > 0
            if isfield(p, 'release_surface') && p.release_surface
                rel = p.k_off * accessible(M);
            else
                rel = p.k_off * M;
            end
        else
            rel = 0;
        end
        m = max(mGuess, 1e-12);
        for it = 1:50
            F = p.D_m * (L * m) + p.lam + rel - a .* m - b .* m.^2;
            if norm(F, inf) < newtTol * max(1, p.lam), break; end
            J  = p.D_m * L - spdiags(a + 2 * b .* m, 0, N, N);
            dm = J \ (-F);
            % damped update with positivity floor
            s = 1;
            for ls = 1:30
                mt = m + s * dm;
                if all(mt > 0), break; end
                s = s / 2;
            end
            m = max(m + s * dm, 1e-12);
        end
        nNewton = nNewton + it;
        maxIter = max(maxIter, it);
    end

% ---- right-hand side -------------------------------------------------------
    function dy = rhs(~, y)
        nRHS = nRHS + 1;
        if qss
            T = y(1:N); P = y(N+1:2*N); M = y(2*N+1:3*N);
            if oligo, S = y(3*N+1:4*N); end
            m = solve_monomer(T, P, M, mWarm);
            mWarm = m;
        else
            T = y(1:N); m = y(N+1:2*N); P = y(2*N+1:3*N); M = y(3*N+1:4*N);
            if oligo, S = y(4*N+1:5*N); end
        end
        A     = accessible(M);
        if hist
            Q = y(end-N+1:end);
            A = Q ./ (1 + (max(M, 0) / p.M_c).^(1/3));
        end
        m2    = m.^2;
        if isfield(p, 'k_off') && p.k_off > 0
            if isfield(p, 'release_surface') && p.release_surface
                netel = 2 * p.k_E * m .* P - p.k_off * A;   % exchange at the shell
            else
                netel = 2 * p.k_E * m .* P - p.k_off * M;   % exchange with all mass
            end
        else
            netel = 2 * p.k_E * m .* P;
        end
        Rtau  = p.k_tau * T .* m2;
        R2    = p.k_2 * A .* m2;
        birth = Rtau + R2;

        dT = p.D_T * (L * T) + p.r_T * T .* (1 - T);
        if oligo
            conv = p.k_conv * S;
            dS = p.D_S * (L * S) + birth - conv;
            dP = conv - p.mu_P * P;
            dM = netel + p.nu * conv - p.mu_M * M;
        else
            dP = birth - p.mu_P * P;
            dM = netel + p.nu * birth - p.mu_M * M;
        end
        if qss
            dy = [dT; dP; dM];
        else
            dm = p.D_m * (L * m) + p.lam * (1 - m) - netel - p.nu * birth;
            dy = [dT; dm; dP; dM];
        end
        if oligo, dy = [dy; dS]; end
        if hist
            sig = p.sig_lo + (1 - p.sig_lo) * max(m, 0);
            dy  = [dy; sig .* max(dM, 0) - p.mu_M * Q];
        end
    end

% ---- integrate -------------------------------------------------------------
if qss
    y0 = [T0; P0; M0];
else
    y0 = [T0; m0; P0; M0];
end
if oligo, y0 = [y0; zeros(N, 1)]; end
if hist,  y0 = [y0; zeros(N, 1)]; end
% T, P and M stay non-negative under these dynamics; A(M) is guarded at M = 0,
% so no NonNegative constraint is imposed (it slows ode15s appreciably here).
odeopt = odeset('RelTol', RelTol, 'AbsTol', AbsTol);
tic;
[t, Y] = ode15s(@rhs, tspan, y0, odeopt);
wall = toc;

% ---- unpack ----------------------------------------------------------------
nt = numel(t);
sol.x = x; sol.t = t(:);
sol.T = Y(:, 1:N);
if qss
    sol.P = Y(:, N+1:2*N);
    sol.M = Y(:, 2*N+1:3*N);
    sol.m = zeros(nt, N);
    mG = ones(N, 1);
    for k = 1:nt
        mk = solve_monomer(sol.T(k,:).', sol.P(k,:).', sol.M(k,:).', mG);
        sol.m(k, :) = mk.';  mG = mk;
    end
else
    sol.m = Y(:, N+1:2*N);
    sol.P = Y(:, 2*N+1:3*N);
    sol.M = Y(:, 3*N+1:4*N);
end

if oligo
    k0 = (3 + ~qss) * N;
    sol.S = Y(:, k0+1:k0+N);
    sol.oligo_length = sqrt(p.D_S / p.k_conv);
end
if hist
    sol.Q = Y(:, end-N+1:end);      % site-weighted mass; capacity uses this
end
sol.p = p; sol.opts = struct('N', N, 'qss', qss, 'RelTol', RelTol, ...
                             'AbsTol', AbsTol, 'newtTol', newtTol);
sol.stats = struct('wall_s', wall, 'nRHS', nRHS, 'nNewton', nNewton, ...
                   'maxNewtonIter', maxIter, 'h', h);
if verbose
    fprintf('[%s] N=%d qss=%d  %.1f s, %d RHS calls, max %d Newton iters\n', ...
            p.name, N, qss, wall, nRHS, maxIter);
end
end

function v = ternary(c, f, g)
if c, v = f(); else, v = g(); end
end
