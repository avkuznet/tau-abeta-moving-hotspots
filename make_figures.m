function make_figures(outdir)
%MAKE_FIGURES  Generate every figure in the manuscript except the schematic.
%
%   make_figures            writes to ../figures
%   make_figures(outdir)    writes to outdir
%
% Produces, at 600 dpi and no wider than 6.5 in, with Arial throughout:
%
%   fig2  space-time maps of T and R_2, reference run, with the R_2 maximum
%   fig3  profiles of R_2 and m at selected times, reference run
%   fig4  verification of the hotspot condition, Eq. (18)
%   fig5  sensitivity to monomer redistribution
%   fig6  the three observable maxima
%   fig7  calibrated run on a 30 mm segment, in physical units
%   fig8  survival of spatial structure against oligomer transport
%   fig9  reversible deposition and history-dependent catalyticity
%
% Fig. 1 is a schematic and is not generated here.

if nargin < 1, outdir = fullfile('..', 'figures'); end
if ~exist(outdir, 'dir'), mkdir(outdir); end
set(groot, 'defaultAxesFontName', 'Arial', 'defaultTextFontName', 'Arial', ...
           'defaultAxesFontSize', 9, 'defaultLineLineWidth', 1.1);
save_fig = @(f, n) exportgraphics(f, fullfile(outdir, [n '.png']), 'Resolution', 600);

% ---------------------------------------------------------------- solutions
fprintf('solving reference case...\n');
p   = hotspot_params('reference');
sol = hotspot_solve(p, struct('N', 161, 'tspan', 0:0.1:30, 'qss', false));
dg  = hotspot_post(sol);
x = sol.x; t = sol.t;

%% Fig 2 --------------------------------------------------------------------
f = figure('Units', 'inches', 'Position', [1 1 6.5 5.2], 'Color', 'w');
ax1 = subplot(2, 1, 1);
imagesc(t, x, sol.T.'); set(gca, 'YDir', 'normal'); colormap(ax1, parula);
c = colorbar; c.Label.String = 'tau catalyst  {\itT}';
ylabel('{\itx} / {\itL}'); title('(a)  effective tau catalyst', 'FontWeight', 'normal');
ax2 = subplot(2, 1, 2);
imagesc(t, x, dg.R2.'); set(gca, 'YDir', 'normal'); colormap(ax2, parula);
c = colorbar; c.Label.String = 'secondary nucleation  {\itR}_2';
hold on; plot(t, dg.x_R2, 'w-', 'LineWidth', 1.6);
plot(t, dg.x_front, 'w--', 'LineWidth', 1.1);
hold off; xlabel('time  {\itt}'); ylabel('{\itx} / {\itL}');
title('(b)  homotypic secondary-nucleation rate', 'FontWeight', 'normal');
legend({'arg max_{\itx} {\itR}_2', 'tau front'}, 'TextColor', 'w', ...
       'Color', 'none', 'EdgeColor', 'w', 'Location', 'southeast');
save_fig(f, 'fig2'); close(f);

%% Fig 3 --------------------------------------------------------------------
tsel = [5 10 15 20 30];
f = figure('Units', 'inches', 'Position', [1 1 6.5 5.0], 'Color', 'w');
cmap = parula(numel(tsel) + 1);
subplot(2, 1, 1); hold on
for i = 1:numel(tsel)
    k = find(t >= tsel(i), 1); plot(x, dg.R2(k, :), 'Color', cmap(i, :));
end
plot(dg.x_R2(ismember(round(t,6), tsel)), ...
     dg.v_R2(ismember(round(t,6), tsel)), 'k.', 'MarkerSize', 9);
hold off; box on; ylabel('{\itR}_2'); title('(a)  secondary-nucleation rate', 'FontWeight', 'normal');
legend(arrayfun(@(v) sprintf('{\\itt} = %g', v), tsel, 'uni', 0), 'Location', 'northwest');
subplot(2, 1, 2); hold on
for i = 1:numel(tsel)
    k = find(t >= tsel(i), 1); plot(x, sol.m(k, :), 'Color', cmap(i, :));
end
hold off; box on; xlabel('{\itx} / {\itL}'); ylabel('{\itm}');
title('(b)  free monomer', 'FontWeight', 'normal'); ylim([0 1.05]);
save_fig(f, 'fig3'); close(f);

%% Fig 4 --------------------------------------------------------------------
tchk = [12 15 20 25 30];
lhs = interp1(t, dg.lhs17, tchk); rhs = interp1(t, dg.rhs17, tchk);
f = figure('Units', 'inches', 'Position', [1 1 4.2 4.0], 'Color', 'w');
lim = [min([lhs rhs]) max([lhs rhs])] + [-1 1] * 0.08 * range([lhs rhs]);
plot(lim, lim, 'k--'); hold on
scatter(rhs, lhs, 46, tchk, 'filled', 'MarkerEdgeColor', 'k');
hold off; box on; axis square; xlim(lim); ylim(lim);
xlabel('-2 \partial_{\itx} ln {\itm}'); ylabel('\partial_{\itx} ln {\itA}');
c = colorbar; c.Label.String = 'time  {\itt}';
title(sprintf('mean departure %.2f%%', 100*mean(abs(lhs-rhs)./abs(rhs))), 'FontWeight', 'normal');
save_fig(f, 'fig4'); close(f);

%% Fig 5 --------------------------------------------------------------------
Dms = [0.01 0.02 0.05 0.10 0.20];
f = figure('Units', 'inches', 'Position', [1 1 6.5 5.0], 'Color', 'w');
cmap = parula(numel(Dms) + 1); lg = cell(1, numel(Dms));
for i = 1:numel(Dms)
    q = hotspot_params('reference'); q.D_m = Dms(i);
    s = hotspot_solve(q, struct('N', 161, 'tspan', 0:0.1:30, 'qss', false));
    g = hotspot_post(s);
    subplot(2, 1, 1); hold on; plot(g.t, g.x_R2, 'Color', cmap(i, :));
    subplot(2, 1, 2); hold on; plot(g.t, g.v_R2 ./ max(g.R2(:, 1), eps), 'Color', cmap(i, :));
    lg{i} = sprintf('{\\itD_m} = %.2f', Dms(i));
end
subplot(2, 1, 1); box on; ylabel('arg max_{\itx} {\itR}_2');
title('(a)  position of the maximum', 'FontWeight', 'normal');
legend(lg, 'Location', 'northwest');
subplot(2, 1, 2); box on; xlabel('time  {\itt}'); ylabel('{\itR}_2 at max / {\itR}_2 at source');
title('(b)  rate at the maximum relative to the source', 'FontWeight', 'normal');
set(gca, 'YScale', 'log');
save_fig(f, 'fig5'); close(f);

%% Fig 6  -------------------------------------------------------------------
f = figure('Units', 'inches', 'Position', [1 1 6.5 3.6], 'Color', 'w');
plot(t, dg.x_birth, '-', 'Color', [0.85 0.33 0.10]); hold on
plot(t, dg.x_R2,   '-', 'Color', [0 0.45 0.74], 'LineWidth', 1.6);
plot(t, dg.x_dMdt, '-', 'Color', [0.47 0.67 0.19], 'LineWidth', 1.6);
plot(t, dg.x_front, 'k--');
hold off; box on; xlabel('time  {\itt}'); ylabel('{\itx} / {\itL}'); ylim([-0.03 1]);
legend({'birth rate {\itR}_\tau + {\itR}_2', 'secondary nucleation {\itR}_2', ...
        'mass accumulation \partial{\itM}/\partial{\itt}', 'tau front'}, ...
       'Location', 'northwest');
title('positions of the three observable maxima, reference run', 'FontWeight', 'normal');
save_fig(f, 'fig6'); close(f);

%% Fig 7  -------------------------------------------------------------------
fprintf('solving calibrated case on a 30 mm segment...\n');
Lc = 30; Nc = round(Lc / 0.05) + 1;
pc  = hotspot_params('calibrated', Lc);
sc  = hotspot_solve(pc, struct('N', Nc, 'tspan', 0:0.25:30, 'qss', true));
gc  = hotspot_post(sc);
mo  = gc.phys.t_months;
f = figure('Units', 'inches', 'Position', [1 1 6.5 5.0], 'Color', 'w');
subplot(2, 1, 1);
plot(mo, gc.phys.x_R2_mm, '-', 'Color', [0 0.45 0.74], 'LineWidth', 1.6); hold on
plot(mo, gc.phys.x_dMdt_mm, '-', 'Color', [0.47 0.67 0.19], 'LineWidth', 1.6);
plot(mo, gc.phys.x_front_mm, 'k--'); hold off; box on
ylabel('distance from source (mm)');
legend({'{\itR}_2 maximum', '\partial{\itM}/\partial{\itt} maximum', 'tau front'}, ...
       'Location', 'northwest');
title('(a)  dimensionally mapped run, 30 mm segment', 'FontWeight', 'normal');
subplot(2, 1, 2);
plot(mo, gc.phys.x_R2_mm - gc.phys.x_front_mm, 'k-'); box on
xlabel('time (months)'); ylabel('lead (mm)');
title('(b)  lead of the {\itR}_2 maximum over the {\itT} = 0.1 contour', 'FontWeight', 'normal');
save_fig(f, 'fig7'); close(f);

%% Fig 8  -------------------------------------------------------------------
% Does the hotspot survive if the nucleation products are mobile? Peter et al.
% (2026) showed that spatial transmission in vitro is carried by diffusing
% oligomers. Scanning the oligomer transport length sqrt(D_S/k_conv) locates
% the value above which spatial structure is lost.
fprintf('scanning oligomer transport length...\n');
MONTH = 720; tauC = 0.50; Lc8 = 30; N8 = 601;
kc = 3.3e-2 * MONTH * tauC;                  % Peter 2026 conversion rate
o8 = struct('N', N8, 'tspan', 0:0.5:30, 'qss', true);
g8 = hotspot_post(hotspot_solve(hotspot_params('calibrated', Lc8), o8));
R0 = g8.R2(end, :); base = max(R0) / R0(1);
Lo = [0.01 0.03 0.05 0.08 0.10 0.13 0.16 0.20 0.25 0.30 0.5 1.0 2.64];
con = zeros(size(Lo));
for i = 1:numel(Lo)
    q = hotspot_params('calibrated', Lc8);
    q.k_conv = kc; q.D_S = kc * (Lo(i) / Lc8)^2;
    gg = hotspot_post(hotspot_solve(q, o8));
    R = gg.R2(end, :); con(i) = max(R) / R(1);
end
f = figure('Units', 'inches', 'Position', [1 1 6.5 3.8], 'Color', 'w');
ax = axes(f); hold(ax, 'on'); yl = [1 2.05];
% shaded: Koffie et al. (2009) halo, 6.5 +/- 4.1 um beyond the plaque core
patch([0.0024 0.0106 0.0106 0.0024], [yl(1) yl(1) yl(2) yl(2)], [0.80 0.90 0.80], ...
      'EdgeColor', 'none', 'FaceAlpha', 0.7);
plot([2.64 2.64], yl, '-', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.4);
plot([0.002 4], [base base], 'k:', 'LineWidth', 1);
plot(Lo, con, 'o-', 'Color', [0 0.45 0.74], 'MarkerFaceColor', 'w', 'LineWidth', 1.3);
hold(ax, 'off'); box on; set(ax, 'XScale', 'log'); ylim(yl); xlim([0.002 4]);
xlabel('oligomer transport length  ({\itD_S}/{\itk}_{conv})^{1/2}  (mm)');
ylabel({'contrast of the maximum', '{\itR}_2 at maximum / {\itR}_2 at source'});
text(0.005, 1.95, {'oligomer halo', '(fixed tissue)'}, 'HorizontalAlignment', 'center', 'FontSize', 8);
text(2.64, 1.95, {'capillary', 'in vitro'}, 'HorizontalAlignment', 'center', ...
     'FontSize', 8, 'Color', [0.85 0.33 0.10]);
title('spatial structure survives only for immobile nucleation products', ...
      'FontWeight', 'normal');
ax.Toolbar = [];
save_fig(f, 'fig8'); close(f);

%% Fig 9  -------------------------------------------------------------------
% The two kinetic refinements the data most clearly demand act in opposite
% directions: reversible deposition suppresses the contrast, history-dependent
% site density raises it.
fprintf('scanning reversibility and history dependence...\n');
Lc9 = 30; N9 = 601; o9 = struct('N', N9, 'tspan', 0:0.5:30, 'qss', true);
sb = hotspot_solve(hotspot_params('calibrated', Lc9), o9);
gb = hotspot_post(sb); Rb = gb.R2(end, :); cb = max(Rb) / Rb(1);
lam9 = sb.p.lam; Msrc = sb.M(end, 1);

relf = [0 0.1 0.25 0.5 0.75 1.0];       % release flux / production flux at the source
koff = relf * lam9 / Msrc;
c1 = zeros(size(koff));
for i = 1:numel(koff)
    q = hotspot_params('calibrated', Lc9);
    if koff(i) > 0, q.k_off = koff(i); end
    R = getfield(hotspot_post(hotspot_solve(q, o9)), 'R2'); R = R(end, :);
    c1(i) = max(R) / R(1);
end
sig = [1 0.5 0.25 0.119 0.05];
c2 = zeros(size(sig)); c3 = zeros(size(sig));
for i = 1:numel(sig)
    q = hotspot_params('calibrated', Lc9);
    if sig(i) < 1, q.sig_lo = sig(i); end
    R = getfield(hotspot_post(hotspot_solve(q, o9)), 'R2'); R = R(end, :);
    c2(i) = max(R) / R(1);
    q.k_off = 0.5 * lam9 / Msrc;        % release at half the production flux
    R = getfield(hotspot_post(hotspot_solve(q, o9)), 'R2'); R = R(end, :);
    c3(i) = max(R) / R(1);
end

f = figure('Units', 'inches', 'Position', [1 1 6.5 3.2], 'Color', 'w');
ax = subplot(1, 2, 1);
plot(relf, c1, 'o-', 'Color', [0.85 0.33 0.10], 'MarkerFaceColor', 'w', 'LineWidth', 1.3);
hold on; plot([1 1], [1 2.7], 'k--'); plot([0 1.1], [1 1], 'k:'); hold off
box on; xlim([0 1.1]); ylim([0.95 2.7]);
xlabel({'release flux / production flux', 'at the tau source'});
ylabel('contrast of the maximum');
title('(a)  reversible deposition', 'FontWeight', 'normal');
text(1.0, 2.5, {'release =', 'production'}, 'HorizontalAlignment', 'right', 'FontSize', 8);
ax.Toolbar = [];

ax = subplot(1, 2, 2);
plot(sig, c2, 'o-', 'Color', [0 0.45 0.74], 'MarkerFaceColor', 'w', 'LineWidth', 1.3); hold on
plot(sig, c3, 's-', 'Color', [0.47 0.67 0.19], 'MarkerFaceColor', 'w', 'LineWidth', 1.3);
plot([0.119 0.119], [0.95 2.7], 'k--'); plot([0 1], [cb cb], 'k:'); hold off
box on; set(ax, 'XDir', 'reverse'); xlim([0.04 1]); ylim([0.95 2.7]);
xlabel({'site density of wake-grown fibrils', '\sigma_{lo} (relative to background)'});
legend({'irreversible', 'with reversibility'}, 'Location', 'northwest', 'FontSize', 8);
title('(b)  history-dependent catalyticity', 'FontWeight', 'normal');
text(0.119, 1.03, 'Hu et al.', 'HorizontalAlignment', 'center', 'FontSize', 8);
ax.Toolbar = [];
save_fig(f, 'fig9'); close(f);

fprintf('figures written to %s\n', outdir);
end
