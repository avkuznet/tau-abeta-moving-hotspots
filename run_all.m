function run_all(mode)
%RUN_ALL  Reproduce the manuscript's numerical results.
%
%   run_all           full report (a few minutes)
%   run_all('quick')  short horizon, coarse grid; for checking the code runs
%
% Sections printed:
%   1  reference case, time-dependent monomer, against the published values
%   2  quasi-steady against fully time-dependent monomer
%   3  the three observable maxima (Section 3.5)
%   4  dimensionally mapped ("calibrated") parameter set (Sections 2.6, 5.6)
%   5  sensitivity to monomer redistribution (Section 5.4)
%   6  grid convergence (Section 4)
%   7  oligomer transport, reversibility and history dependence (5.7, 5.8)
%
% Conventions, as in the manuscript: positions quoted at "grid resolution"
% are the grid point of largest R_2; hotspot_post returns parabolic sub-grid
% positions. Speeds are least-squares slopes over t = 10-30 (reference) and
% 5-15 months (mapped runs).

if nargin < 1, mode = 'full'; end
quick = strcmpi(mode, 'quick');
if quick
    Nref = 81; tsp = 0:0.5:12;
else
    Nref = 161; tsp = 0:0.1:30;
end
tag = @(s) fprintf('\n%s\n%s\n', s, repmat('-', 1, numel(s)));

%% 1 ------------------------------------------------------------------------
tag('1. REFERENCE CASE, TIME-DEPENDENT MONOMER, AGAINST THE PUBLISHED VALUES');
fprintf(['  The quasi-steady reduction is NOT used here. At the reference\n' ...
         '  parameters lambda/r_T = 1, so monomer and tau evolve on the same\n' ...
         '  scale and Eq. (16) does not apply; see Section 2.\n']);
p   = hotspot_params('reference');
sol = hotspot_solve(p, struct('N', Nref, 'tspan', tsp, 'qss', false, 'verbose', true));
dg  = hotspot_post(sol);
xg  = grid_argmax(dg);

pub_t = [10 12 15 20 25 30];
pub_x = [0.394 0.525 0.663 0.756 0.794 0.813];
fprintf('\n  R_2 maximum position\n');
fprintf('  %6s %10s %10s %10s %10s\n', 't', 'grid', 'sub-grid', 'published', 'months');
for k = 1:numel(pub_t)
    if pub_t(k) > max(sol.t), continue; end
    j = nearest(dg.t, pub_t(k));
    fprintf('  %6.1f %10.4f %10.4f %10.3f %10.2f\n', pub_t(k), xg(j), dg.x_R2(j), pub_x(k), ...
            pub_t(k) * p.scale.tau_months);
end
jd = find(xg > 0, 1);
fprintf('\n  first output time with the maximum off x = 0   %6.2f   (published 8.3)\n', dg.t(jd));
fprintf('  same, interpolated sub-grid (x_R2 > 2h)        %6.2f\n', dg.t_detach);
fprintf('  Da_EC at the source at that time               %6.2f   (published 15.6)\n', ...
        2 * p.k_E * sol.P(jd, 1) / p.lam);
fprintf('  source R_2 peak at t                           %6.2f   (published 3.1)\n', dg.t_peak_source);
kb = find(dg.n_peaks >= 2, 1);
fprintf('  R_2 first bimodal at t                         %6.2f\n', dg.t(kb));
ko = find(dg.x_front > dg.x_R2 & dg.t > dg.t(jd), 1);
if ~isempty(ko), fprintf('  T = 0.1 contour overtakes the maximum at t     %6.2f\n', dg.t(ko)); end
fprintf('  hotspot speed, LS over t = 10-30               %7.4f  (published 0.0171)\n', ...
        lsqw(dg.t, dg.x_R2, 10, 30));
fprintf('  Fisher front speed 2 sqrt(D_T r_T)             %7.4f  (published 0.0548)\n', dg.v_front_theory);
report_at(dg, xg, 12, p);
[eg, es] = err18(dg, sol, xg, 10, 30);
fprintf('  Eq. (18) relative error, t = 10-30: grid maximum mean %.2f%%, max %.2f%%;\n', ...
        100 * mean(eg), 100 * max(eg));
fprintf('                                      sub-grid maximum mean %.3f%%, max %.3f%%\n', ...
        100 * mean(es), 100 * max(es));
pubDa = [2.50 15.0; 7.55 45.3];
for q = 1:2
    tt = 10 + 20 * (q - 1); if tt > max(dg.t), continue; end
    j = nearest(dg.t, tt); i = find(dg.x == xg(j), 1);
    sink = 2 * p.k_E * sol.P(j, i);
    fprintf(['  t = %2d: Da_EC = %.2f, Da_D = %.1f at the grid maximum (published %.2f, %.1f);' ...
             ' %.2f, %.1f at the sub-grid maximum\n'], tt, sink / p.lam, sink / p.D_m, ...
            pubDa(q, 1), pubDa(q, 2), dg.Da_EC(j), dg.Da_D(j));
end
fprintf('  in months on the 7.2 mm reference segment: maximum leaves the source at %.2f\n', ...
        dg.t(jd) * p.scale.tau_months);

%% 2 ------------------------------------------------------------------------
tag('2. QUASI-STEADY AGAINST TIME-DEPENDENT MONOMER');
solF = hotspot_solve(p, struct('N', Nref, 'tspan', tsp, 'qss', true, 'verbose', true));
dgF  = hotspot_post(solF);
fprintf('\n  %6s %12s %12s %10s\n', 't', 'x_R2 (full)', 'x_R2 (qss)', 'diff');
for tt = pub_t
    if tt > max(sol.t), continue; end
    a = interp1(dg.t, dg.x_R2, tt); b = interp1(dgF.t, dgF.x_R2, tt);
    fprintf('  %6.1f %12.4f %12.4f %10.4f\n', tt, a, b, a - b);
end
fprintf('  max |m_qss - m_full| over all times and nodes: %.3e   (published 1.4e-1)\n', ...
        max(abs(sol.m(:) - solF.m(:))));

%% 3 ------------------------------------------------------------------------
tag('3. THE THREE OBSERVABLE MAXIMA  (Section 3.5)');
fprintf('  %6s %8s %10s %10s %11s %12s %6s\n', 't', 'months', 'x(R_2)', 'x(dM/dt)', ...
        'x(birth)', 'x(tau front)', 'peaks');
for tt = pub_t
    if tt > max(sol.t), continue; end
    fprintf('  %6.1f %8.2f %10.3f %10.3f %11.3f %12.3f %6d\n', tt, tt*p.scale.tau_months, ...
            interp1(dg.t, dg.x_R2, tt), interp1(dg.t, dg.x_dMdt, tt), ...
            interp1(dg.t, dg.x_birth, tt), interp1(dg.t, dg.x_front, tt), ...
            round(interp1(dg.t, dg.n_peaks, tt)));
end
fprintf('  max x(dM/dt) over the run: %.4f   (the mass-accumulation maximum stays at the source)\n', ...
        max(dg.x_dMdt));

%% 4 ------------------------------------------------------------------------
tag('4. DIMENSIONALLY MAPPED PARAMETER SET  (Sections 2.6, 5.6)');
% A fixed PHYSICAL resolution, h = 0.05 mm, resolves the 0.66 mm tau front;
% the initial tau pulse is 0.288 mm wide on every segment (hotspot_params).
Lcal = 30; if quick, Lcal = 20; end
[solC, dgC] = mapped_run(Lcal, tsp, true);
mo = dgC.phys.t_months;
fprintf('\n  %d mm segment\n', Lcal);
fprintf('  %8s %10s %10s %14s %14s %12s\n', 'months', 'x_R2 (mm)', 'x_dMdt', 'Da_EC at max', ...
        'Da_EC source', 'm at source');
DaS = 2 * solC.p.k_E * solC.P(:, 1) / solC.p.lam;
for mm = [1 2 4 5 6 8 10 12.5 15]
    if mm > max(mo), continue; end
    j = nearest(mo, mm);
    fprintf('  %8.1f %10.2f %10.2f %14.2f %14.2f %12.3f\n', mm, dgC.phys.x_R2_mm(j), ...
            dgC.phys.x_dMdt_mm(j), dgC.Da_EC(j), DaS(j), dgC.m_at_source(j));
end
jd = find(grid_argmax(dgC) > 0, 1);
fprintf('\n  maximum leaves the source at   %7.2f months\n', mo(jd));
fprintf('  source R_2 peak at             %7.2f months\n', dgC.phys.t_peak_source_mo);
fprintf('  hotspot speed, LS 5-15 months  %7.3f mm/month\n', lsqw(mo, dgC.phys.x_R2_mm, 5, 15));
fprintf('  T = 0.1 contour speed          %7.3f mm/month   (2 sqrt(D_T r_T))\n', dgC.phys.v_front_mm_mo);
fprintf('  lead over T = 0.1 at %.0f months %6.2f mm\n', mo(end), ...
        dgC.phys.x_R2_mm(end) - dgC.phys.x_front_mm(end));
fprintf('  monomer at source at %.0f months %6.3f of background (%.0f%% depletion)\n', mo(end), ...
        dgC.m_at_source(end), 100 * (1 - dgC.m_at_source(end)));
fprintf('  max x(dM/dt) over the run      %7.3f mm\n', max(dgC.phys.x_dMdt_mm));
R = dgC.R2(end, :);
fprintf('  contrast R_2(max)/R_2(source)  %7.3f at %.0f months\n', max(R) / R(1), mo(end));
if ~quick
    fprintf('\n  segment length   speed LS 5-15 mo   x_R2 at 15 mo   lead at 15 mo\n');
    for Lx = [20 30 40 50]
        [~, gx] = mapped_run(Lx, tsp, false);
        fprintf('  %8d mm %14.3f %15.2f %15.2f\n', Lx, lsqw(gx.phys.t_months, gx.phys.x_R2_mm, 5, 15), ...
                gx.phys.x_R2_mm(end), gx.phys.x_R2_mm(end) - gx.phys.x_front_mm(end));
    end
end
fprintf('\n  quasi-steady check at the mapped set (lambda/r_T = 422):\n');
solCF = hotspot_solve(solC.p, struct('N', numel(solC.x), 'tspan', tsp, 'qss', false));
fprintf('    max |m_qss - m_full| = %.2e   (reference set: %.2e)\n', ...
        max(abs(solC.m(:) - solCF.m(:))), max(abs(sol.m(:) - solF.m(:))));

%% 5 ------------------------------------------------------------------------
tag('5. SENSITIVITY TO MONOMER REDISTRIBUTION');
Dms = [0.01 0.02 0.05 0.10 0.20];
pubD = [4.4 5.5 8.3 12.8 NaN];
pubX = [0.900 0.858 0.758 0.625 NaN];
fprintf('  %6s %16s %12s %16s %12s\n', 'D_m', 'detach (grid)', 'published', 'x_h(20) (grid)', 'published');
for i = 1:numel(Dms)
    q = hotspot_params('reference'); q.D_m = Dms(i);
    s = hotspot_solve(q, struct('N', Nref, 'tspan', tsp, 'qss', false));
    g = hotspot_post(s); gx = grid_argmax(g);
    j = find(gx > 0, 1); td = NaN; if ~isempty(j), td = g.t(j); end
    x20 = NaN; if max(s.t) >= 20, x20 = gx(nearest(g.t, 20)); end
    fprintf('  %6.2f %16.2f %12.1f %16.3f %12.3f\n', Dms(i), td, pubD(i), x20, pubX(i));
end

%% 6 ------------------------------------------------------------------------
tag('6. GRID CONVERGENCE');
Ns = [81 161 321];
store = cell(numel(Ns), 1); grd = cell(numel(Ns), 1);
for i = 1:numel(Ns)
    s = hotspot_solve(p, struct('N', Ns(i), 'tspan', tsp, 'qss', false));
    store{i} = hotspot_post(s); grd{i} = grid_argmax(store{i});
    fprintf('  N = %3d done (%.1f s)\n', Ns(i), s.stats.wall_s);
end
chk = pub_t(pub_t <= max(tsp));
fprintf('\n  %6s %12s %12s %12s\n', 't', 'N=81', 'N=161', 'N=321');
for tt = chk
    fprintf('  %6.1f', tt);
    for i = 1:numel(Ns), fprintf(' %12.4f', interp1(store{i}.t, store{i}.x_R2, tt)); end
    fprintf('\n');
end
dg_ = 0; ds_ = 0; dp_ = 0;
for tt = chk
    j = nearest(store{2}.t, tt);
    dg_ = max(dg_, abs(grd{2}(j) - grd{3}(j)));
    ds_ = max(ds_, abs(store{2}.x_R2(j) - store{3}.x_R2(j)));
    dp_ = max(dp_, abs(max(store{2}.R2(j, :)) - max(store{3}.R2(j, :))) / max(store{3}.R2(j, :)));
end
fprintf('\n  161 -> 321: max change in x_h %.4f at grid resolution (published 0.0031),\n', dg_);
fprintf('              %.4f sub-grid (published 0.0006), peak R_2 %.3f%% (published < 0.1%%)\n', ...
        ds_, 100 * dp_);

%% 7 ------------------------------------------------------------------------
if ~quick
    tag('7. OLIGOMER TRANSPORT, REVERSIBILITY, HISTORY DEPENDENCE  (5.7, 5.8)');
    sections_57_58();
end
fprintf('\nDone.\n');
end

% ---------------------------------------------------------------------------
function [s, g] = mapped_run(Lmm, tsp, verbose)
pc = hotspot_params('calibrated', Lmm);
s  = hotspot_solve(pc, struct('N', round(Lmm / 0.05) + 1, 'tspan', tsp, 'qss', true, 'verbose', verbose));
g  = hotspot_post(s);
end

function sections_57_58()
% Same settings as make_figures (Figs. 8 and 9): 30 mm, N = 601, t = 0:0.5:30.
Lc = 30; o = struct('N', 601, 'tspan', 0:0.5:30, 'qss', true);
MONTH = 720; kc = 3.3e-2 * MONTH * 0.5;          % Peter 2026 conversion rate, per time unit
sb = hotspot_solve(hotspot_params('calibrated', Lc), o); gb = hotspot_post(sb);
c0 = contrast(gb); x0 = gb.x_R2(end) * Lc;
fprintf('  base model: contrast %.3f, maximum at %.2f mm (15 months)\n', c0, x0);
fprintf('\n  oligomer transport length (mm)   contrast   maximum (mm)\n');
for Lo = [0.10 0.13 0.15 0.16 0.20 0.25 0.30 0.50 1.00 2.64]
    q = hotspot_params('calibrated', Lc); q.k_conv = kc; q.D_S = kc * (Lo / Lc)^2;
    g = hotspot_post(hotspot_solve(q, o));
    fprintf('  %20.2f %14.3f %13.2f\n', Lo, contrast(g), g.x_R2(end) * Lc);
end
q = hotspot_params('calibrated', Lc); q.k_conv = kc; q.D_S = 1e-14;
g = hotspot_post(hotspot_solve(q, o));
fprintf('  conversion delay only (D_S -> 0): contrast %.3f, maximum %.2f mm (shift %.1f%%)\n', ...
        contrast(g), g.x_R2(end) * Lc, 100 * abs(g.x_R2(end) * Lc - x0) / x0);
lam = sb.p.lam; Ms = sb.M(end, 1); As = gb.A(end, 1);
fprintf('\n  release flux / production flux at the source    contrast   maximum (mm)\n');
xs = [];
for rf = [0.5 1.0]
    q = hotspot_params('calibrated', Lc); q.k_off = rf * lam / Ms;
    g = hotspot_post(hotspot_solve(q, o)); xs(end+1) = g.x_R2(end) * Lc; %#ok<AGROW>
    fprintf('  %8.1f, whole-mass exchange %24.3f %13.2f\n', rf, contrast(g), xs(end));
    q.k_off = rf * lam / As; q.release_surface = true;
    g = hotspot_post(hotspot_solve(q, o));
    fprintf('  %8.1f, surface-limited exchange %20.3f %13.2f\n', rf, contrast(g), g.x_R2(end) * Lc);
end
fprintf('\n  sigma_lo   contrast (irreversible)   contrast (release = production/2)\n');
for sg = [0.5 0.25 0.119 0.05]
    q = hotspot_params('calibrated', Lc); q.sig_lo = sg;
    g1 = hotspot_post(hotspot_solve(q, o));
    q.k_off = 0.5 * lam / Ms;
    g2 = hotspot_post(hotspot_solve(q, o));
    xs(end+1:end+2) = [g1.x_R2(end), g2.x_R2(end)] * Lc;
    fprintf('  %8.3f %18.3f %26.3f\n', sg, contrast(g1), contrast(g2));
end
fprintf('  largest shift of the maximum across these combinations: %.1f%%\n', ...
        100 * max(abs(xs - x0)) / x0);
fprintf('  (sigma_lo = 0.119 is 99/834, the Brichos-site reduction of Hu et al. 2026)\n');
end

function c = contrast(g)
R = g.R2(end, :); c = max(R) / R(1);
end

function xg = grid_argmax(dg)
[~, i] = max(dg.R2, [], 2); xg = dg.x(i); xg = xg(:);
end

function j = nearest(t, tt)
[~, j] = min(abs(t - tt));
end

function s = lsqw(t, y, t1, t2)
w = t >= t1 & t <= t2; t = t(w); y = y(w); t = t(:); y = y(:);
s = ((t - mean(t)).' * (y - mean(y))) / sum((t - mean(t)).^2);
end

function [eg, es] = err18(dg, sol, xg, t1, t2)
% relative residual of Eq. (18) at the grid maximum and at the sub-grid maximum
h = dg.x(2) - dg.x(1); kk = find(dg.t >= t1 & dg.t <= t2);
eg = zeros(numel(kk), 1); es = dg.err17(kk);
lnA = log(max(dg.A, realmin)); lnm = log(max(sol.m, realmin));
for q = 1:numel(kk)
    k = kk(q); i = find(dg.x == xg(k), 1);
    dA = gradient(lnA(k, :), h); dm = gradient(lnm(k, :), h);
    eg(q) = abs(dA(i) + 2 * dm(i)) / abs(2 * dm(i));
end
end

function report_at(dg, xg, tt, p)
if tt > max(dg.t), return; end
j = nearest(dg.t, tt); i = find(dg.x == xg(j), 1);
fprintf('\n  at t = %g (%.1f months), grid maximum x = %.4f\n', tt, tt * p.scale.tau_months, xg(j));
m_hot = sqrt(dg.R2(j, i) / (p.k_2 * dg.A(j, i)));     % m at the grid maximum
fprintf('    m  source/hot  %7.3f / %7.3f   (published 0.105 / 0.386)\n', dg.m_at_source(j), m_hot);
fprintf('    A  source/hot  %7.3f / %7.3f   (published 1.976 / 0.283)\n', dg.A_at_source(j), dg.A(j, i));
fprintf('    R_2 hot/source %7.3f\n', dg.R2(j, i) / dg.R2(j, 1));
end
