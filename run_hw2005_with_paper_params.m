%% RUN_HW2005_WITH_PAPER_PARAMS.M
%  Replicates Table II of the paper using the Table III structural parameters.
%  Pipeline: build state space -> VFI -> simulate firm panel -> compute
%  moments -> print comparison table
%
%  Expected runtime:
%    'test'  mode: 1-2 minutes (smoke test)
%    'paper' mode: 1-3 hours (full replication of Table II)
%
%  Reference: Hennessy & Whited (2005), Journal of Finance, Vol.60, No.3
%
%  Changes relative to the original version:
%    1. Loop variable renamed from sim to sim_val, to avoid shadowing the
%       MATLAB built-in function sim
%    2. Added a panel sanity check before computing moments, to catch
%       NaN/Inf early

clear; clc; close all;
rng(42);   % fixed random seed for reproducibility

fprintf('============================================================\n');
fprintf('  Hennessy & Whited (2005) "Debt Dynamics" Replication\n');
fprintf('  Objective: reproduce Table II simulated moments using Table III parameters\n');
fprintf('============================================================\n\n');

%% ================================================================
%  1. Select run mode
%     'test'  - minimal grid, quick sanity check
%     'paper' - paper-standard grid, full replication
%% ================================================================
RUN_MODE = 'paper';   % switch mode by editing this line

[fixed, vfi_opts, sim_opts, ~] = config_hw2005(RUN_MODE);
fprintf('\n');

%% ================================================================
%  2. Table III structural parameters
%
%  Parameter ordering (consistent across all sub-functions):
%    params = [alpha, s, delta, rho, sigma_eps, mu_tau, sigma_tau, lambda]
%    index:       1     2    3     4      5         6        7         8
%% ================================================================
%                 alpha    s     delta   rho   sigma_eps  mu_tau  sigma_tau  lambda
b_paper = [      0.551,  0.592,  0.100,  0.740,  0.123,  -12.267,  9.246,   0.059];

param_names = {'alpha','s','delta','rho','sigma_eps','mu_tau','sigma_tau','lambda'};
param_desc  = {'Curvature of profit function', 'Fire-sale discount', 'Capital depreciation rate', ...
               'Persistence of shock', 'Std. dev. of shock', 'Location of tax rate', ...
               'Scale of tax rate', 'Equity issuance cost'};

fprintf('Table III structural parameters:\n');
fprintf('  %-12s  %9s  %s\n', 'Parameter', 'Value', 'Description');
fprintf('  %s\n', repmat('-', 50, 1));
for pp = 1:8
    fprintf('  %-12s  %9.4f  %s\n', param_names{pp}, b_paper(pp), param_desc{pp});
end
fprintf('\n');

%% ================================================================
%  3. Data moments M_hat (Table II "Actual Moments")
%     Order matches the output of compute_moments.m exactly
%% ================================================================
M_hat = [
    0.079;   % Moment 1:  Mean(I/A)
    0.006;   % Moment 2:  Var(I/A)
    0.146;   % Moment 3:  Mean(EBITDA/A)
    0.075;   % Moment 4:  Mean(net debt/A)
    0.042;   % Moment 5:  Mean(equity issuance/A)
    0.099;   % Moment 6:  Freq(equity issuance)
    0.019;   % Moment 7:  Investment-Q sensitivity
   -0.080;   % Moment 8:  Debt-Q sensitivity
    0.583;   % Moment 9:  Serial correlation of income/assets
    0.117;   % Moment 10: SD of shock to income/assets
];

% Weighting matrix (diagonal, inverse of squared data moments)
% Note: this is a validation weight only, not the paper's influence-function
% weights, so the objective value is not directly comparable to the
% paper's chi2 = 4.906.
W = diag(1 ./ (M_hat.^2 + 1e-8));

%% ================================================================
%  4. Validate the corporate_tax function
%% ================================================================
fprintf('--- Validating corporate_tax ---\n');
[tc0,  g0]  = corporate_tax(0,     b_paper(6), b_paper(7));
[tc_hi, ~]  = corporate_tax(1e4,   b_paper(6), b_paper(7));
[tc_lo, ~]  = corporate_tax(-1e4,  b_paper(6), b_paper(7));
expected_tc0 = 0.35 * normcdf(-b_paper(6)/b_paper(7));
fprintf('  tau_c(0)    = %.4f  (theoretical = %.4f)\n', tc0, expected_tc0);
fprintf('  tau_c(+inf) = %.4f  (should approach 0.3500)\n', tc_hi);
fprintf('  tau_c(-inf) = %.4f  (should approach 0.0000)\n', tc_lo);
fprintf('  g(0)        = %.2e  (must equal 0)\n', g0);
assert(abs(g0) < 1e-12, 'corporate_tax: g(0) != 0, function is incorrect');
fprintf('  corporate_tax validation passed\n\n');

%% ================================================================
%  5. Build the state space
%% ================================================================
fprintf('--- Step 1: Building the state space ---\n');
t_grid = tic;
grids = build_state_space(b_paper, fixed);
fprintf('  elapsed: %.1f s\n\n', toc(t_grid));

%% ================================================================
%  6. Value Function Iteration
%% ================================================================
fprintf('--- Step 2: Value function iteration (VFI) ---\n');
fprintf('  state space: %d x %d x %d = %d states\n', ...
        fixed.Nk, fixed.Np, fixed.Nz, fixed.Nk*fixed.Np*fixed.Nz);
fprintf('  convergence criterion: tol=%.0e, max_iter=%d\n', vfi_opts.tol, vfi_opts.max_iter);
t_vfi = tic;
[V, k_idx_pol, p_idx_pol, vfi_history] = solve_vfi_fast(grids, b_paper, fixed, vfi_opts);
fprintf('  VFI total time: %.1f s\n\n', toc(t_vfi));

% VFI sanity checks
assert(~any(isinf(V(:))),  'VFI: V contains Inf, check solve_vfi_fast');
assert(~any(isnan(V(:))),  'VFI: V contains NaN, check solve_vfi_fast');
assert(min(V(:)) > -1e10,  'VFI: V has an implausibly small minimum (< -1e10), possibly not converged');

%% ================================================================
%  7. Simulate the firm panel
%% ================================================================
fprintf('--- Step 3: Simulating the firm panel ---\n');
fprintf('  setup: S=%d panels x N=%d firms x T=%d periods, keep last %d periods\n', ...
        sim_opts.S_sim, sim_opts.N_firm, sim_opts.T_tot, sim_opts.T_keep);
t_sim = tic;
panel = simulate_firms(V, k_idx_pol, p_idx_pol, grids, b_paper, fixed, sim_opts);
fprintf('  simulation total time: %.1f s\n\n', toc(t_sim));

%% ================================================================
%  7.5 Panel sanity check (before computing moments, to catch problems early)
%% ================================================================
fprintf('--- Panel sanity check ---\n');
fields_to_check = {'I_A','EBITDA_A','Q','debt_A','eq_A','income_A'};
all_ok = true;
for ff = 1:numel(fields_to_check)
    fn  = fields_to_check{ff};
    vec = panel.(fn);
    n_nan = sum(isnan(vec));
    n_inf = sum(isinf(vec));
    if n_nan > 0 || n_inf > 0
        fprintf('  [WARNING] panel.%-12s: %d NaN, %d Inf\n', fn, n_nan, n_inf);
        all_ok = false;
    end
end
% Range plausibility checks
checks = {
    'I/A in (-1,5)',      panel.I_A,      -1,   5;
    'EBITDA/A in (0,2)',  panel.EBITDA_A,  0,   2;
    'Q in (0,100)',        panel.Q,         0, 100;
    'debt/A in (-2,2)',   panel.debt_A,   -2,   2;
    'eq_A >= 0',          panel.eq_A,      0, Inf;
};
for cc = 1:size(checks,1)
    label = checks{cc,1};
    vec   = checks{cc,2};
    lo    = checks{cc,3};
    hi    = checks{cc,4};
    n_out = sum(vec < lo | vec > hi);
    if n_out > 0
        pct_out = 100*n_out/numel(vec);
        fprintf('  [WARNING] %-20s: %d out-of-range (%.1f%%)\n', label, n_out, pct_out);
        if pct_out > 5, all_ok = false; end
    end
end
if all_ok
    fprintf('  panel checks passed\n');
end
fprintf('\n');

%% ================================================================
%  8. Compute simulated moments
%% ================================================================
fprintf('--- Step 4: Computing simulated moments ---\n');
m_sim = compute_moments(panel);
fprintf('  done\n\n');

% Moment sanity checks
assert(~any(isnan(m_sim)), 'compute_moments returned NaN, check panel data');
assert(numel(m_sim) == 10, 'compute_moments did not return 10 moments');

%% ================================================================
%  9. Objective function value
%% ================================================================
G   = M_hat - m_sim;
obj = G' * W * G;

%% ================================================================
%  10. Print result table (corresponding to Table II)
%% ================================================================
moment_names = {
    'Mean(I/A)';
    'Var(I/A)';
    'Mean(EBITDA/A)';
    'Mean(net debt/A)';
    'Mean(eq issuance/A)';
    'Freq(eq issuance)';
    'Investment-Q sensitivity';
    'Debt-Q sensitivity';
    'Serial corr(income/A)';
    'SD(shock income/A)'
};

fprintf('\n');
fprintf('%s\n', repmat('=', 72, 1));
fprintf('  TABLE II REPLICATION  [mode: %s]\n', upper(RUN_MODE));
fprintf('%s\n', repmat('=', 72, 1));
fprintf('%-26s  %10s  %10s  %10s  %8s\n', ...
        'Moment', 'Actual', 'Simulated', 'Diff', 'Diff%');
fprintf('%s\n', repmat('-', 72, 1));

max_diff_pct = 0;
for mm = 1:10
    act     = M_hat(mm);
    sim_val = m_sim(mm);          % renamed: sim -> sim_val, to avoid shadowing a built-in
    dif     = sim_val - act;
    if abs(act) > 1e-8
        pct = 100 * dif / abs(act);
    else
        pct = NaN;
    end
    if ~isnan(pct)
        max_diff_pct = max(max_diff_pct, abs(pct));
    end

    flag = '';
    if ~isnan(pct) && abs(pct) > 20, flag = ' **'; end
    if ~isnan(pct) && abs(pct) > 5 && abs(pct) <= 20, flag = ' *'; end

    fprintf('%-26s  %10.4f  %10.4f  %10.4f  %7.1f%%%s\n', ...
            moment_names{mm}, act, sim_val, dif, pct, flag);
end
fprintf('%s\n', repmat('-', 72, 1));
fprintf('  GMM objective value  = %.4f\n', obj);
fprintf('  max relative deviation = %.1f%%\n', max_diff_pct);
fprintf('  Note: diagonal weighting matrix; not directly comparable to the paper''s chi2 = 4.906\n');
fprintf('%s\n', repmat('=', 72, 1));

%% ================================================================
%  11. Additional diagnostics: basic statistics of simulated data
%% ================================================================
fprintf('\n--- Simulated data diagnostics ---\n');
fprintf('  total observations:       %8d\n',   length(panel.I_A));
fprintf('  equity issuance freq:      %8.3f (paper: 0.099)\n', mean(panel.eq_dummy));
fprintf('  mean net debt/assets:      %8.3f (paper: 0.075)\n', mean(panel.debt_A));
fprintf('  fraction of saving firms:  %8.3f (paper approx: 0.10)\n', mean(panel.debt_A < 0));
fprintf('  mean Tobin Q:              %8.3f (paper approx: 1-2)\n',  mean(panel.Q));
fprintf('  mean I/A:                  %8.3f (paper: 0.079)\n',  mean(panel.I_A));
fprintf('  mean EBITDA/A:             %8.3f (paper: 0.146)\n',  mean(panel.EBITDA_A));

%% ================================================================
%  12. Save results
%% ================================================================
save_name = sprintf('hw2005_%s_results.mat', RUN_MODE);
save(save_name, 'b_paper', 'M_hat', 'm_sim', 'obj', 'panel', ...
     'grids', 'V', 'k_idx_pol', 'p_idx_pol', '-v7.3');
fprintf('\nResults saved to: %s\n', save_name);

%% ================================================================
%  13. Plot simulated moments against data moments
%% ================================================================
try
    figure('Name', 'Table II Moments Comparison', 'Position', [100, 100, 900, 500]);

    x = 1:10;
    bar_data = [M_hat, m_sim];
    b_h = bar(x, bar_data, 0.7);
    b_h(1).FaceColor = [0.2, 0.4, 0.8];
    b_h(2).FaceColor = [0.9, 0.4, 0.2];

    set(gca, 'XTick', x, 'XTickLabel', ...
        {'I/A','Var(I/A)','EBITDA/A','Debt/A','EqIss/A',...
         'FreqEq','IQ-sens','DQ-sens','AR(1)','SD(shock)'}, ...
        'XTickLabelRotation', 30, 'FontSize', 9);
    legend({'Actual (Table II)', 'Simulated'}, 'Location', 'northeast');
    title(sprintf('HW2005 Table II: Actual vs Simulated Moments [%s mode]', RUN_MODE));
    ylabel('Moment Value');
    grid on; box on;

    saveas(gcf, sprintf('hw2005_%s_moments.png', RUN_MODE));
    fprintf('Comparison figure saved to: hw2005_%s_moments.png\n', RUN_MODE);
catch
    fprintf('  (figure output skipped)\n');
end

fprintf('\nAll done.\n');
