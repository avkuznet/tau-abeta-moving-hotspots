function p = hotspot_params(name, Lmm)
%HOTSPOT_PARAMS  Parameter sets for the moving-hotspot reaction-diffusion model.
%
%   p = HOTSPOT_PARAMS('reference')   Table 1 of the manuscript. Reproduces the
%                                     published proof-of-principle results.
%   p = HOTSPOT_PARAMS('calibrated')  Monomer sector set from measurement; see
%                                     Section 2.6 and the notes below.
%
% Every field is dimensionless. The mapping to physical units is fixed by the
% two tau-front parameters and is returned in p.scale:
%
%   r_T   = 0.30 (model)  <->  0.6 month^-1      (Meisl et al. 2021, P301S mouse)
%   D_T   = 0.0025        <->  ~1e-13 m^2 s^-1   (Meisl et al. 2021)
%       =>  time unit  tau = 0.50 months
%           length unit L   = 7.20 mm
%
% so the reference horizon t = 0..30 is 15 months across 7.2 mm of tissue.
%
% CALIBRATED SET
%   lam  and D_m are set from measurement:
%       lambda_m = ln2/1.97 h = 9.77e-5 s^-1  (Cirrito et al. 2003)  -> 126.7
%       D_m      = 62.3 um^2 s^-1             (Waters 2010)          ->   1.557
%   k_E is then chosen so that the kinetic Damkohler number reached in the
%   depleted wake matches the measured value, Da_EC = 2k_E P/lambda = 1.4 to
%   2.3 (Hong et al. 2011, from tracer recovery and from steady-state ISF
%   Abeta). Because the absolute scale of the fibril-number variable P is
%   arbitrary in the dimensionless model, that scale is absorbed into k_E;
%   k_tau, k_2, mu_P, nu and M_c are left at their reference values, no
%   independent measurement constraining them.
%
%   k_E = 15 was found by sweeping k_E and reading off Da_EC. The calibration
%   is weakly determined: depletion at the source varies only from 87% to 97%
%   as k_E runs from 10 to 509, because a smaller elongation sink leaves more
%   monomer, which raises nucleation and hence P, largely cancelling the
%   change. Da_EC is the more sensitive target and is what k_E is set by.
%
%   Achieved at k_E = 15, for checking the run:
%       Da_EC               1.0 to 2.5 over the run   (measured 1.4 to 2.3)
%       screening length    500 to 790 um             (estimated 442 to 530)
%       depletion at source 81.1% at 6 months         (Hong 2011, Abeta42: 81%)
%       monomer relaxation is ~1.4e3 times faster than tau replication, which
%       is what justifies the quasi-steady reduction of Eq. (16).

if nargin < 1, name = 'reference'; end

switch lower(name)
    case 'reference'
        p.D_T   = 0.0025;   % tau spreading coefficient
        p.r_T   = 0.30;     % tau amplification rate
        p.lam   = 0.30;     % monomer restoration/clearance
        p.D_m   = 0.05;     % monomer redistribution
        p.k_E   = 4.0;      % elongation
        p.k_tau = 0.50;     % heterotypic tau-Abeta nucleation
        p.k_2   = 0.80;     % homotypic secondary nucleation
        p.mu_P  = 0.02;     % fibril-number loss
        p.mu_M  = 0.0;      % fibril-mass loss
        p.nu    = 0.02;     % mass per nucleation event
        p.M_c   = 0.20;     % shielding crossover mass

    case 'calibrated'
        % The time unit is fixed by r_T alone: 0.30 / (0.6 month^-1) = 0.50
        % months. The length unit is then free, and is taken as an argument so
        % that the domain can be made long enough for the maximum to stay
        % interior. Lmm = 7.20 recovers the value implied by Table 1's D_T.
        % D_T and D_m are rescaled with L so that the physical spreading
        % coefficient, monomer diffusivity, front speed and screening length
        % are unchanged; lam, k_E and the nucleation coefficients are rates
        % and so do not depend on L.
        if nargin < 2 || isempty(Lmm), Lmm = 7.20; end
        tau_months = 0.50;
        D_T_phys   = 0.2592;   % mm^2/month, from 1e-13 m^2/s   (Meisl 2021)
        D_m_phys   = 161.5;    % mm^2/month, from 62.3 um^2/s   (Waters 2010)

        p.r_T   = 0.30;
        p.D_T   = D_T_phys * tau_months / Lmm^2;
        p.lam   = 126.67;   % lambda_m * tau, Cirrito 2003
        p.D_m   = D_m_phys * tau_months / Lmm^2;
        p.k_E   = 15.0;     % set by calibration, see note below
        p.k_tau = 0.50;
        p.k_2   = 0.80;
        p.mu_P  = 0.02;
        p.mu_M  = 0.0;
        p.nu    = 0.02;
        p.M_c   = 0.20;

    otherwise
        error('hotspot_params:unknown', 'Unknown parameter set "%s".', name);
end

p.name = lower(name);

% initial tau pulse, T(x,0) = T0*exp(-(x/w)^2)
p.T0 = 0.02;
p.w  = 0.04;

% physical scales (Section 2.6)
p.scale.tau_months = 0.50;   % one time unit, months
if strcmp(p.name,'calibrated'), p.scale.L_mm = Lmm; else, p.scale.L_mm = 7.20; end
p.scale.note = 'tau and L fixed by r_T and D_T against Meisl et al. (2021)';
end
