function [fixed, vfi_opts, sim_opts, sa_opts] = config_hw2005(mode)
% CONFIG_HW2005  Configuration for the Hennessy & Whited (2005) replication.
%
%   Key design principle:
%     fixed contains only parameters that are truly fixed across modes
%     (tax rates, interest rate, grid sizes).
%     The 8 structural parameters of Table III are passed separately to
%     each sub-function through the b_paper vector, and are deliberately
%     NOT placed in fixed. Otherwise, during estimation b would vary while
%     fixed stays constant, silently corrupting the objective.
%
%   Mode options:
%     'test'  - minimal grid, 1-2 minutes, smoke test only
%     'fast'  - coarse grid, 15-30 minutes, results indicative only
%     'paper' - paper-standard grid, 1-3 hours, required to replicate Table II
%     'full'  - same grid as paper, used for SMM estimation
%
%   OUTPUTS
%     fixed    - struct: r, tau_d, tau_i, Nk, Np, Nz
%     vfi_opts - struct: tol, max_iter
%     sim_opts - struct: S_sim, N_firm, T_tot, T_keep
%     sa_opts  - struct: simulated-annealing options (full mode only)

if nargin < 1, mode = 'paper'; end

%% ----------------------------------------------------------------
%  Truly fixed parameters (shared by all modes; from Sections I and V
%  of the paper)
%% ----------------------------------------------------------------
fixed.r     = 0.025;   % real risk-free rate (paper, Section V)
fixed.tau_d = 0.12;    % dividend tax rate (Graham 2000; paper, Section V)
fixed.tau_i = 0.25;    % personal tax on interest income (slightly below Graham's 0.296; paper, Section V)

% NOTE: alpha, s, delta, rho, sigma_eps, mu_tau, sigma_tau, lambda
% are the structural parameters to be estimated; they are passed through
% the params/b vector and are intentionally NOT defined here.

%% ----------------------------------------------------------------
%  Grid sizes and solver options by mode
%% ----------------------------------------------------------------
switch lower(mode)

    %% ---- Test mode (minimal, 1-2 minutes) ----
    case 'test'
        % Smallest feasible grid, used only to verify the code runs
        fixed.Nk = 7;
        fixed.Np = 12;
        fixed.Nz = 7;

        vfi_opts.tol      = 1e-3;
        vfi_opts.max_iter = 300;

        sim_opts.S_sim  = 1;
        sim_opts.N_firm = 200;
        sim_opts.T_tot  = 20;
        sim_opts.T_keep = 9;

        sa_opts = struct();   % no estimation in test mode

    %% ---- Fast mode (coarse grid, 15-30 minutes) ----
    case 'fast'
        % Deviates from the paper grid; indicative results only,
        % not for replicating Table II
        fixed.Nk = 15;
        fixed.Np = 30;
        fixed.Nz = 15;

        vfi_opts.tol      = 1e-5;
        vfi_opts.max_iter = 800;

        sim_opts.S_sim  = 6;
        sim_opts.N_firm = 10000;
        sim_opts.T_tot  = 50;
        sim_opts.T_keep = 9;

        sa_opts.MaxIterations          = 30;
        sa_opts.MaxFunctionEvaluations = 200;
        sa_opts.InitialTemperature     = 10;

    %% ---- Paper-standard mode (required to replicate Table II) ----
    case 'paper'
        % The paper specifies Nz = 20 and Nk = 21 (j = 0,...,20);
        % Np is not specified in the paper, 100 points are sufficient
        fixed.Nk = 21;
        fixed.Np = 100;
        fixed.Nz = 20;

        vfi_opts.tol      = 1e-6;
        vfi_opts.max_iter = 2000;

        % Paper, Section VI: S = 6 panels, N = 10000 firms, T = 50 periods,
        % last 9 periods retained
        sim_opts.S_sim  = 6;
        sim_opts.N_firm = 10000;
        sim_opts.T_tot  = 50;
        sim_opts.T_keep = 9;

        sa_opts = struct();   % no estimation in validation mode

    %% ---- Full estimation mode (SMM) ----
    case 'full'
        % Same grid as paper mode, plus optimization options
        fixed.Nk = 21;
        fixed.Np = 100;
        fixed.Nz = 20;

        vfi_opts.tol      = 1e-6;
        vfi_opts.max_iter = 2000;

        sim_opts.S_sim  = 6;
        sim_opts.N_firm = 10000;
        sim_opts.T_tot  = 50;
        sim_opts.T_keep = 9;

        sa_opts.MaxIterations          = 5000;
        sa_opts.MaxFunctionEvaluations = 50000;
        sa_opts.InitialTemperature     = 100;
        sa_opts.ReannealingInterval    = 100;

    otherwise
        error('config_hw2005: mode must be test/fast/paper/full, got "%s"', mode);
end

%% ----------------------------------------------------------------
%  Print configuration summary
%% ----------------------------------------------------------------
fprintf('[config] mode="%s" | Nk=%d Np=%d Nz=%d | tol=%.0e max_iter=%d\n', ...
        mode, fixed.Nk, fixed.Np, fixed.Nz, vfi_opts.tol, vfi_opts.max_iter);
fprintf('[config] simulation: S=%d N=%d T=%d keep=%d | r=%.3f tau_d=%.2f tau_i=%.2f\n', ...
        sim_opts.S_sim, sim_opts.N_firm, sim_opts.T_tot, sim_opts.T_keep, ...
        fixed.r, fixed.tau_d, fixed.tau_i);

end
