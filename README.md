# Debt Dynamics — Replication of Hennessy & Whited (2005)

A MATLAB implementation of **Hennessy and Whited (2005), "Debt Dynamics,"** *Journal of Finance*, Vol. 60, No. 3, pp. 1129–1165. The code solves the dynamic model of the levered firm with endogenous investment, financing, savings, and equity issuance, simulates a panel of firms, and computes the corresponding simulated moments, following the model and the structural parameters (Table III) of the paper.

## Files

All MATLAB files below live in the `code/` folder. Pre-run outputs (results and figures) are in `results/`.

| File | Purpose |
|---|---|
| `code/run_hw2005_with_paper_params.m` | Main script. Builds the state space, runs value function iteration (VFI), simulates the firm panel, computes moments, prints a Table II comparison, saves results, and plots actual vs. simulated moments. |
| `code/config_hw2005.m` | Configuration for four modes: `test`, `fast`, `paper`, `full`. |
| `code/build_state_space.m` | Constructs the capital, debt, and productivity grids, including the collateral upper bound. |
| `code/solve_vfi_fast.m` | Value function iteration with policy-function updates; returns convergence history. |
| `code/simulate_firms.m` | Simulates S panels of N firms over T periods and records the panel used for moments. |
| `code/compute_moments.m` | Computes the 10 simulated moments (investment, debt, equity issuance, Q sensitivities, AR(1) and shock SD of income). |
| `code/corporate_tax.m` | Convex marginal tax rate and tax bill as specified in the paper (Assumption 6). |
| `code/tauchen86.m` | Discretization of the AR(1) productivity process (Tauchen 1986). |
| `code/diagnose_Q.m` | Optional diagnostic script for the investment–Q sensitivity. |
| `code/generate_report_figures.m` | Generates report figures: VFI convergence, policy distributions, summary table, policy functions, and the Table II moment comparison. |

## Requirements

- MATLAB R2021b or later (the code uses `bsxfun`, `yline`, `sgtitle`, and `-v7.3` MAT-file saving).
- No third-party toolboxes required (Statistics and Machine Learning Toolbox functions such as `normcdf`/`normpdf` are used by `corporate_tax.m`).

## How to run

Open MATLAB, `cd` into the `code/` folder, and run:

```matlab
run_hw2005_with_paper_params
```

The main script is configured with `RUN_MODE = 'paper'`, which replicates the paper's grid and is the mode reported below.

### Run modes

| Mode | Grid (Nk × Np × Nz) | Runtime | Use |
|---|---|---|---|
| `test` | 7 × 12 × 7 | ~1–2 min | Smoke test |
| `fast` | 15 × 30 × 15 | ~15–30 min | Quick, indicative |
| `paper` | 21 × 100 × 20 | ~1–3 h | Full replication of Table II |
| `full` | 21 × 100 × 20 | long | SMM estimation (simulated annealing) |

To switch modes, edit `RUN_MODE` in `code/run_hw2005_with_paper_params.m`.

After the main script finishes, report figures can be generated with:

```matlab
generate_report_figures
```

This script loads the variables from the workspace (`V`, `k_idx_pol`, `p_idx_pol`, `vfi_history`, `grids`, `panel`, `m_sim`, `M_hat`, `b_paper`, `fixed`), so it must be run in the same session, or after loading the saved `.mat` results.

## Results

The paper-standard grid (`paper` mode, `rng(42)`) solves the model and produces the full set of simulated moments, saved to `results/hw2005_paper_results.mat` together with the policy functions, value function, grids, and panel. A comparison table is printed at the end of the run, and the moments are shown against the paper's Table II in `results/hw2005_paper_moments.png`.

**Replication fidelity.** The model solution and simulation pipeline are stable and reproducible (identical moments across runs with a fixed seed), but the simulated moments deviate from the moments reported in Table II of the paper. The largest gaps are in debt and equity issuance (`Mean(net debt/A)` ≈ 0.010 vs. 0.127 in the paper; `Freq(eq issuance)` ≈ 0.019 vs. 0.091) and in the income dynamics (`SD(shock income/A)` ≈ 0.027 vs. 0.117). These gaps are discussed in the next section and largely reflect deliberate departures in how the moments are constructed on the short simulated panel. This repository should therefore be read as a faithful implementation of the model's numerical machinery, not as a claim to have reproduced the paper's exact moment values.

**Notes on comparability with the paper:**

- The diagonal weighting matrix `W = diag(1 ./ (M_hat.^2 + 1e-8))` is used for validation only. It is not the influence-function weighting used in the paper, so the GMM objective value is not directly comparable to the paper's reported `chi2 = 4.906`.

## Implementation notes and departures from a strict reading of the paper

The following choices were made during development and are deliberately kept, as they improve numerical behavior on the short simulated panel:

1. **Moments 9–10 (serial correlation and shock SD):** a fixed-effects demeaned direct OLS AR(1) regression is used instead of a differenced IV estimator. On a short panel (T = 9), the differenced IV estimator showed severe small-sample bias (negative estimates despite clear positive autocorrelation in the data). The shock SD is estimated directly from the AR(1) residuals as `std(u_it)`.
2. **Q calculation in the simulation:** follows the paper's definition strictly, `Q = (EV + b'/(1+r)) / k`, using the expected continuation value under the next-period policy and the *current-period* capital in the denominator (fixes a sign reversal found during diagnostics).
3. **Tauchen grid width:** `m_std = 2` instead of 3 in `tauchen86.m`. For a persistent process (`rho = 0.74`) with N = 20 points, `m = 3` makes the endpoint probabilities tiny and crowds the interior points; `m = 2` covers the 95% unconditional distribution much better.
4. The loop variable in the main script is named `sim_val` rather than `sim`, to avoid shadowing MATLAB's built-in `sim`.

## Reproducibility

- Random seed is fixed (`rng(42)`) at the start of the main script.
- No randomization inside the VFI; only the firm-panel simulation draws shocks.
- Results are saved as `results/hw2005_paper_results.mat` (MAT-file v7.3).

## Reference

Hennessy, C. A., & Whited, T. M. (2005). Debt Dynamics. *The Journal of Finance*, 60(3), 1129–1165. https://doi.org/10.1111/j.1540-6261.2005.00762.x

## License

For personal and educational use. If you use this code in your own work, please cite the original paper.
