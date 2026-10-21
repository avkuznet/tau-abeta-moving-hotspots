function run_all(mode)
%RUN_ALL  Reproduce the manuscript's numerical results.
%
%   run_all           full report
%   run_all('quick')  short horizon, coarse grid; for checking the code runs
%
% Sections printed:
%   1  reference case, time-dependent monomer, against the published values
%   2  quasi-steady against fully time-dependent monomer
%   3  the three observable maxima (Section 3.5)
%   4  calibrated parameter set in physical units (Section 2.6)
%   5  sensitivity to monomer redistribution
%   6  grid convergence

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

pub_t = [10 12 15 20 25 30];
pub_x = [0.394 0.525 0.663 0.756 0.794 0.813];
fprintf('\n  R_2 maximum position\n');
fprintf('  %6s %10s %10s %9s %12s\n', 't', 'computed', 'published', 'diff', 'months');
for k = 1:numel(pub_t)
    if pub_t(k) > max(sol.t), continue; end
    xc = interp1(dg.t, dg.x_R2, pub_t(k));
    fprintf('  %6.1f %10.3f %10.3f %9.3f %12.2f\n', pub_t(k), xc, pub_x(k), ...
            xc - pub_x(k), pub_t(k) * p.scale.tau_months);
end
fprintf('\n  detachment time        %7.2f   (published 8.3)\n', dg.t_detach);
fprintf('  source R_2 peak at t   %7.2f   (published 3.1)\n', dg.t_peak_source);
fprintf('  hotspot speed          %7.4f   (published 0.0171)\n', dg.v_hotspot);
fprintf('  Fisher front speed     %7.4f   (published 0.0548)\n', dg.v_front_theory);
report_at(dg, 12, p);
fprintf('  Eq. (18) relative error over the moving window: max %.3f%%, mean %.3f%%\n', ...
        100*max(dg.err17(isfinite(dg.err17) & dg.x_R2 > 0.05)), ...
        100*mean(dg.err17(isfinite(dg.err17) & dg.x_R2 > 0.05)));

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
fprintf('  max |m_qss - m_full| over all times and nodes: %.3e\n', ...
        max(abs(sol.m(:) - solF.m(:))));
fprintf('  wall clock: qss %.1f s, full %.1f s\n', sol.stats.wall_s, solF.stats.wall_s);

%% 3 ------------------------------------------------------------------------
tag('3. THE THREE OBSERVABLE MAXIMA  (Section 3.5)');
fprintf('  R_2 is not measurable in tissue; dM/dt is what PET integrates and\n');
fprintf('  R_tau+R_2 is the rate at which new deposits appear.\n\n');
fprintf('  %6s %8s %10s %10s %11s %12s %6s\n', 't', 'months', 'x(R_2)', 'x(dM/dt)', ...
        'x(birth)', 'x(tau front)', 'peaks');
for tt = pub_t
    if tt > max(sol.t), continue; end
    fprintf('  %6.1f %8.2f %10.3f %10.3f %11.3f %12.3f %6d\n', tt, tt*p.scale.tau_months, ...
            interp1(dg.t, dg.x_R2, tt), interp1(dg.t, dg.x_dMdt, tt), ...
            interp1(dg.t, dg.x_birth, tt), interp1(dg.t, dg.x_front, tt), ...
            round(interp1(dg.t, dg.n_peaks, tt)));
end
fprintf('\n  Two results worth noting. The mass-accumulation maximum stays at the\n');
fprintf('  tau source throughout, because dM/dt = 2 k_E m P peaks where P is\n');
fprintf('  largest despite the depletion there, so the quantity amyloid imaging\n');
fprintf('  reports and the quantity the model organises sit at opposite ends of\n');
fprintf('  the domain. And R_2 is bimodal: a peak at the source and a peak ahead\n');
fprintf('  of it exchange dominance, so the maximum does not migrate smoothly.\n');
fprintf('  The birth-rate maximum tracks the R_2 maximum closely, which makes\n');
fprintf('  the rate of appearance of new deposits the usable experimental proxy.\n');
fprintf('\n  in mm (L = %.2f mm):\n', p.scale.L_mm);
fprintf('  %8s %10s %10s %11s %12s\n', 'months', 'x(R_2)', 'x(dM/dt)', 'x(birth)', 'lag of R_2');
for tt = pub_t
    if tt > max(sol.t), continue; end
    fprintf('  %8.2f %10.2f %10.2f %11.2f %12.2f\n', tt*p.scale.tau_months, ...
            interp1(dg.t, dg.phys.x_R2_mm, tt), interp1(dg.t, dg.phys.x_dMdt_mm, tt), ...
            interp1(dg.t, dg.phys.x_birth_mm, tt), interp1(dg.t, dg.phys.lag_mm, tt));
end

%% 4 ------------------------------------------------------------------------
tag('4. CALIBRATED PARAMETER SET  (Section 2.6)');
% The calibrated case needs a domain long enough for the maximum to stay
% interior, and a fixed PHYSICAL resolution: the tau front is 0.66 mm wide,
% so h is held at 0.05 mm rather than N held fixed. L = 30 mm and L = 40 mm
% give the same hotspot speed to three decimals; L = 20 mm does not.
Lcal = 30; Ncal = round(Lcal/0.05) + 1;
if quick, Lcal = 20; Ncal = 201; end
pc   = hotspot_params('calibrated', Lcal);
solC = hotspot_solve(pc, struct('N', Ncal, 'tspan', tsp, 'qss', true, 'verbose', true));
dgC  = hotspot_post(solC);
fprintf('\n  target Da_EC ~ 2.3 (Hong 2011); screening length ~530 um\n');
fprintf('  %8s %10s %10s %10s %12s %10s\n', 'months', 'x_R2 (mm)', 'x_dMdt', 'Da_EC', 'Da_D', 'ell (um)');
for tt = pub_t
    if tt > max(solC.t), continue; end
    fprintf('  %8.2f %10.2f %10.2f %10.2f %12.1f %10.0f\n', tt*pc.scale.tau_months, ...
            interp1(dgC.t, dgC.phys.x_R2_mm, tt), interp1(dgC.t, dgC.phys.x_dMdt_mm, tt), ...
            interp1(dgC.t, dgC.Da_EC, tt), interp1(dgC.t, dgC.Da_D, tt), ...
            interp1(dgC.t, dgC.phys.ell_um, tt));
end
fprintf('\n  detachment            %7.2f months\n', dgC.phys.t_detach_months);
fprintf('  hotspot speed         %7.3f mm/month\n', dgC.phys.v_hotspot_mm_mo);
fprintf('  tau front speed       %7.3f mm/month   (He 2018: 0.6 to 1.2)\n', dgC.phys.v_front_mm_mo);
k81 = find(dgC.m_at_source <= 0.19, 1, 'first');
if ~isempty(k81)
    fprintf('  81%% depletion at      %7.1f months  (Hong 2011, Abeta42: 81%% by 24 mo)\n', ...
            dgC.t(k81) * pc.scale.tau_months);
end
fprintf('  monomer at source     %7.3f of background at end of run\n', dgC.m_at_source(end));

fprintf('\n  quasi-steady check at the calibrated set (lambda/r_T = 422):\n');
solCF = hotspot_solve(pc, struct('N', Ncal, 'tspan', tsp, 'qss', false));
fprintf('    max |m_qss - m_full| = %.2e   (reference set: %.2e)\n', ...
        max(abs(solC.m(:) - solCF.m(:))), max(abs(sol.m(:) - solF.m(:))));

%% 5 ------------------------------------------------------------------------
tag('5. SENSITIVITY TO MONOMER REDISTRIBUTION');
Dms = [0.01 0.02 0.05 0.10 0.20];
pubD = [4.4 5.5 8.3 12.8 NaN];
pubX = [0.900 0.858 0.758 0.625 NaN];
fprintf('  %6s %12s %12s %12s %12s\n', 'D_m', 'detach', 'published', 'x_h(t=20)', 'published');
for i = 1:numel(Dms)
    q = hotspot_params('reference'); q.D_m = Dms(i);
    s = hotspot_solve(q, struct('N', Nref, 'tspan', tsp, 'qss', false));
    g = hotspot_post(s);
    x20 = NaN; if max(s.t) >= 20, x20 = interp1(g.t, g.x_R2, 20); end
    fprintf('  %6.2f %12.2f %12.1f %12.3f %12.3f\n', Dms(i), g.t_detach, pubD(i), x20, pubX(i));
end

%% 6 ------------------------------------------------------------------------
tag('6. GRID CONVERGENCE');
Ns = [81 161 321];
store = cell(numel(Ns), 1);
for i = 1:numel(Ns)
    s = hotspot_solve(p, struct('N', Ns(i), 'tspan', tsp, 'qss', false));
    store{i} = hotspot_post(s);
    fprintf('  N = %3d done (%.1f s)\n', Ns(i), s.stats.wall_s);
end
chk = pub_t(pub_t <= max(tsp));
fprintf('\n  %6s %12s %12s %12s\n', 't', 'N=81', 'N=161', 'N=321');
for tt = chk
    fprintf('  %6.1f', tt);
    for i = 1:numel(Ns), fprintf(' %12.4f', interp1(store{i}.t, store{i}.x_R2, tt)); end
    fprintf('\n');
end
d12 = max(abs(arrayfun(@(tt) interp1(store{2}.t, store{2}.x_R2, tt) - ...
                             interp1(store{3}.t, store{3}.x_R2, tt), chk)));
fprintf('\n  max |x_h(161) - x_h(321)| = %.4f   (published 0.0031)\n', d12);

fprintf('\nDone.\n');
end

function report_at(dg, tt, p)
if tt > max(dg.t), return; end
fprintf('\n  at t = %g (%.1f months)\n', tt, tt * p.scale.tau_months);
fprintf('    x_h            %7.3f\n', interp1(dg.t, dg.x_R2, tt));
fprintf('    m  source/hot  %7.3f / %7.3f   (published 0.105 / 0.386)\n', ...
        interp1(dg.t, dg.m_at_source, tt), interp1(dg.t, dg.m_at_hotspot, tt));
fprintf('    A  source/hot  %7.3f / %7.3f   (published 1.976 / 0.283)\n', ...
        interp1(dg.t, dg.A_at_source, tt), interp1(dg.t, dg.A_at_hotspot, tt));
fprintf('    Da_EC / Da_D   %7.2f / %7.1f\n', ...
        interp1(dg.t, dg.Da_EC, tt), interp1(dg.t, dg.Da_D, tt));
end
