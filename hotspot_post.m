function dg = hotspot_post(sol)
%HOTSPOT_POST  Diagnostics for a solution returned by HOTSPOT_SOLVE.
%
%   dg = HOTSPOT_POST(sol)
%
% Returns, for every output time, the position of three maxima:
%
%   dg.x_R2    argmax_x R_2            the quantity the model organises, but
%                                      not measurable in tissue
%   dg.x_dMdt  argmax_x dM/dt          what amyloid PET integrates regionally
%                                      and what two-photon imaging resolves
%   dg.x_birth argmax_x (R_tau + R_2)  the rate at which new deposits appear
%
% Section 3.5 of the manuscript: these three need not coincide, and the
% distance between them is itself a prediction. Peak positions are refined to
% sub-grid accuracy by fitting a parabola through the three samples straddling
% the discrete maximum.
%
% Also returned: the kinetic and spatial Damkohler numbers at each maximum, the
% monomer screening length, the position of the tau front, a check of the
% hotspot condition Eq. (18), and physical-unit versions of all positions and
% speeds using the calibration in sol.p.scale.

p = sol.p; x = sol.x; t = sol.t;
N = numel(x); nt = numel(t); h = x(2) - x(1);

% Exposed catalytic capacity. With a history-dependent site density the
% solver carries a site-weighted mass Q and the capacity is built from it, so
% the diagnostics must use the same quantity the solve used.
if isfield(sol, 'Q')
    A = sol.Q ./ (1 + (max(sol.M, 0) / p.M_c).^(1/3));
else
    A = max(sol.M, 0) ./ (1 + (max(sol.M, 0) / p.M_c).^(1/3));
end
m2    = sol.m.^2;
Rtau  = p.k_tau * sol.T .* m2;
R2    = p.k_2 * A .* m2;
birth = Rtau + R2;
dMdt  = 2 * p.k_E * sol.m .* sol.P + p.nu * birth - p.mu_M * sol.M;

dg.t = t; dg.x = x;
dg.A = A; dg.R2 = R2; dg.Rtau = Rtau; dg.birth = birth; dg.dMdt = dMdt;

[dg.x_R2,    dg.v_R2]    = peaks(x, R2);
[dg.x_dMdt,  dg.v_dMdt]  = peaks(x, dMdt);
[dg.x_birth, dg.v_birth] = peaks(x, birth);

% tau front: outermost crossing of T = thr. A Fisher front is only established
% once T saturates near the source, which with T(x,0) = 0.02 and r_T = 0.30
% takes t ~ ln(1/0.02)/r_T ~ 13; a level set well below 1/2 is therefore
% tracked so that the leading edge is defined over most of the horizon.
if isfield(sol, 'frontLevel'), thr = sol.frontLevel; else, thr = 0.10; end
dg.frontLevel = thr;
dg.T_max = max(sol.T, [], 2);
dg.x_front = nan(nt, 1);
for k = 1:nt
    Tk = sol.T(k, :).';
    j  = find(Tk >= thr, 1, 'last');
    if ~isempty(j) && j < N
        dg.x_front(k) = x(j) + h * (Tk(j) - thr) / (Tk(j) - Tk(j+1));
    elseif ~isempty(j)
        dg.x_front(k) = x(end);
    end
end

% Is the R_2 maximum a single peak that migrates, or a second peak that grows
% ahead of the first and overtakes it? Count interior local maxima above 1% of
% the global maximum and record the two largest.
dg.n_peaks = zeros(nt, 1);
dg.x_peak1 = nan(nt, 1); dg.x_peak2 = nan(nt, 1);
dg.v_peak1 = nan(nt, 1); dg.v_peak2 = nan(nt, 1);
for k = 1:nt
    f = R2(k, :).';
    loc = find([false; (f(2:end-1) > f(1:end-2)) & (f(2:end-1) >= f(3:end)); false]);
    if f(1) > f(2), loc = [1; loc]; end %#ok<AGROW>
    loc = loc(f(loc) > 0.01 * max(f));
    dg.n_peaks(k) = numel(loc);
    [~, ord] = sort(f(loc), 'descend');
    if ~isempty(ord)
        dg.x_peak1(k) = x(loc(ord(1))); dg.v_peak1(k) = f(loc(ord(1)));
    end
    if numel(ord) > 1
        dg.x_peak2(k) = x(loc(ord(2))); dg.v_peak2(k) = f(loc(ord(2)));
    end
end

% Damkohler numbers and screening length, evaluated at the R_2 maximum
dg.Da_EC = nan(nt, 1); dg.Da_D = nan(nt, 1); dg.ell = nan(nt, 1);
dg.m_at_hotspot = nan(nt, 1); dg.A_at_hotspot = nan(nt, 1);
dg.m_at_source  = sol.m(:, 1); dg.A_at_source = A(:, 1);
dg.lhs17 = nan(nt, 1); dg.rhs17 = nan(nt, 1);

for k = 1:nt
    xi = dg.x_R2(k);
    Pk = interp1(x, sol.P(k, :).', xi, 'linear', 'extrap');
    sink = 2 * p.k_E * Pk;
    dg.Da_EC(k) = sink / p.lam;
    dg.Da_D(k)  = sink / p.D_m;              % L = 1 (domain scale)
    dg.ell(k)   = sqrt(p.D_m / max(sink, eps));
    dg.m_at_hotspot(k) = interp1(x, sol.m(k, :).', xi, 'linear', 'extrap');
    dg.A_at_hotspot(k) = interp1(x, A(k, :).',     xi, 'linear', 'extrap');

    % Eq. (18):  d(ln A)/dx  =  -2 d(ln m)/dx  at the interior maximum
    lnA = log(max(A(k, :).', realmin));
    lnm = log(max(sol.m(k, :).', realmin));
    dlnA = gradient(lnA, h);
    dlnm = gradient(lnm, h);
    dg.lhs17(k) = interp1(x, dlnA, xi, 'linear', 'extrap');
    dg.rhs17(k) = -2 * interp1(x, dlnm, xi, 'linear', 'extrap');
end
dg.err17 = abs(dg.lhs17 - dg.rhs17) ./ max(abs(dg.rhs17), eps);

% detachment: first time the R_2 maximum is clear of the source
thr = 2 * h;
kd = find(dg.x_R2 > thr, 1, 'first');
if isempty(kd)
    dg.t_detach = NaN;
else
    if kd > 1   % linear interpolation in t onto x_R2 = thr
        t1 = t(kd-1); t2 = t(kd); a1 = dg.x_R2(kd-1); a2 = dg.x_R2(kd);
        dg.t_detach = t1 + (thr - a1) * (t2 - t1) / (a2 - a1);
    else
        dg.t_detach = t(kd);
    end
end

% source-region R_2 maximum in time
[~, ks] = max(R2(:, 1));
dg.t_peak_source = t(ks);

% speeds by least squares over the window in which the maximum is interior
w = dg.x_R2 > thr & dg.x_R2 < 1 - thr;
dg.v_hotspot = lsq_slope(t(w), dg.x_R2(w));
wf = isfinite(dg.x_front) & dg.x_front < 1 - thr;
dg.v_front_measured = lsq_slope(t(wf), dg.x_front(wf));
dg.v_front_theory   = 2 * sqrt(p.D_T * p.r_T);

% ---- physical units --------------------------------------------------------
Lmm = p.scale.L_mm; tm = p.scale.tau_months;
dg.phys.t_months     = t * tm;
dg.phys.x_R2_mm      = dg.x_R2 * Lmm;
dg.phys.x_dMdt_mm    = dg.x_dMdt * Lmm;
dg.phys.x_birth_mm   = dg.x_birth * Lmm;
dg.phys.x_front_mm   = dg.x_front * Lmm;
dg.phys.ell_um       = dg.ell * Lmm * 1000;
dg.phys.t_detach_months   = dg.t_detach * tm;
dg.phys.t_peak_source_mo  = dg.t_peak_source * tm;
dg.phys.v_hotspot_mm_mo   = dg.v_hotspot * Lmm / tm;
dg.phys.v_front_mm_mo     = dg.v_front_theory * Lmm / tm;
dg.phys.lag_mm            = dg.phys.x_front_mm - dg.phys.x_R2_mm;
end

% ---------------------------------------------------------------------------
function [xp, vp] = peaks(x, F)
% sub-grid peak position and value for every row of F
nt = size(F, 1); xp = zeros(nt, 1); vp = zeros(nt, 1);
h = x(2) - x(1);
for k = 1:nt
    f = F(k, :);
    [vmax, i] = max(f);
    if i == 1 || i == numel(f) || ~isfinite(vmax)
        xp(k) = x(i); vp(k) = vmax; continue;
    end
    y1 = f(i-1); y2 = f(i); y3 = f(i+1);
    den = y1 - 2*y2 + y3;
    if den == 0
        xp(k) = x(i); vp(k) = y2; continue;
    end
    dlt = 0.5 * (y1 - y3) / den;
    if ~isfinite(dlt) || abs(dlt) > 1, dlt = 0; end
    xp(k) = x(i) + dlt * h;
    vp(k) = y2 - 0.25 * (y1 - y3) * dlt;
end
end

function s = lsq_slope(t, y)
if numel(t) < 2, s = NaN; return; end
t = t(:); y = y(:);
s = ((t - mean(t)).' * (y - mean(y))) / sum((t - mean(t)).^2);
end
