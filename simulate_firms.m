function panel = simulate_firms(V, k_idx_pol, p_idx_pol, grids, params, fixed, sim_opts)
% SIMULATE_FIRMS  Firm-panel simulation for HW2005 - with a corrected Q
% calculation (expected continuation value and current capital; fixes the
% sign reversal).

if nargin < 7 || isempty(sim_opts)
    sim_opts.S_sim   = 6;
    sim_opts.N_firm  = 10000;
    sim_opts.T_tot   = 50;
    sim_opts.T_keep  = 9;
end
S_sim   = sim_opts.S_sim;
N_firm  = sim_opts.N_firm;
T_tot   = sim_opts.T_tot;
T_keep  = sim_opts.T_keep;

alpha     = params(1);
s_val     = params(2);
delta     = params(3);
mu_tau    = params(6);
sigma_tau = params(7);
lambda    = params(8);
r         = fixed.r;
tau_d     = fixed.tau_d;

k_grid = grids.k_grid;
p_grid = grids.p_grid;
z_grid = grids.z_grid;
Pi     = grids.Pi;
Nk     = grids.Nk;
Np     = grids.Np;
Nz     = grids.Nz;
Pi_cum = cumsum(Pi, 2);

inv1r = 1 / (1 + r);

total_obs    = S_sim * N_firm * T_keep;
I_A_out      = zeros(total_obs, 1);
EBITDA_A_out = zeros(total_obs, 1);
Q_out        = zeros(total_obs, 1);
debt_A_out   = zeros(total_obs, 1);
eq_A_out     = zeros(total_obs, 1);
eq_dum_out   = false(total_obs, 1);
income_A_out = zeros(total_obs, 1);
firm_id_out  = zeros(total_obs, 1);
year_id_out  = zeros(total_obs, 1);
obs_ptr      = 0;

for sim = 1:S_sim
    fprintf('  panel %d/%d\n', sim, S_sim);
    ik_all = repmat(ceil(Nk/2), N_firm, 1);
    ip_all = repmat(ceil(Np/2), N_firm, 1);
    iz_all = repmat(ceil(Nz/2), N_firm, 1);

    for t = 1:T_tot
        keep = (t > T_tot - T_keep);
        lin_cur = sub2ind([Nk, Np, Nz], ik_all, ip_all, iz_all);
        ikp_all = double(k_idx_pol(lin_cur));
        ipp_all = double(p_idx_pol(lin_cur));

        if keep
            k_cur  = k_grid(ik_all);
            b_cur  = p_grid(ip_all);
            z_cur  = z_grid(iz_all);
            k_next = k_grid(ikp_all);
            b_next = p_grid(ipp_all);

            bv_cur  = b_cur  * inv1r;
            bv_next = b_next * inv1r;

            pi_all = z_cur .* k_cur .^ alpha;
            y_all = pi_all - delta * k_cur - r * bv_cur;
            [~, g_all] = corporate_tax(y_all, mu_tau, sigma_tau);
            pi_net = pi_all - g_all;

            cash_short = (pi_net < bv_cur);
            phi_s = double(cash_short);
            phi_n = 1 - phi_s;
            CF_all = pi_net - bv_cur;
            n_all = phi_s .* max(0, (bv_cur - pi_net) / s_val);
            k_eff = k_cur .* (1 - delta) - n_all;
            invest_all = k_next - k_eff;
            phi_all = phi_n + s_val .* phi_s;
            bracket_all = phi_all .* CF_all - invest_all + bv_next;

            % Keep V_cur for later use in firm value
            V_cur = V(lin_cur);

            % ================== Corrected Q calculation ==================
            % 1. Compute the conditional expected value EV
            V_reshape = reshape(V, [Nk*Np, Nz]);
            EV_matrix = V_reshape * Pi'; 

            % 2. Index the expected value for the next-period policy
            lin_next_state = (ikp_all - 1) * Np + ipp_all; 
            EV_next = zeros(N_firm, 1);
            for f = 1:N_firm
                EV_next(f) = EV_matrix(lin_next_state(f), iz_all(f));
            end

            % 3. Follow the paper definition strictly: the denominator is
            %    the current-period capital k_cur
            Q_all = (EV_next + bv_next) ./ k_cur;
            % =============================================================

            pos = (bracket_all >= 0);
            neg = ~pos;
            eq_iss = zeros(N_firm, 1);
            eq_iss(neg) = -bracket_all(neg);

            firm_val = V_cur + bv_next;
            debt_r = zeros(N_firm, 1);
            ok = (firm_val > 0);
            debt_r(ok) = bv_next(ok) ./ firm_val(ok);
            debt_r(~ok) = bv_next(~ok) ./ k_next(~ok);

            idx = (obs_ptr+1):(obs_ptr+N_firm);
            I_A_out(idx)      = invest_all ./ k_cur;
            EBITDA_A_out(idx) = pi_all     ./ k_cur;
            Q_out(idx)        = Q_all;
            debt_A_out(idx)   = debt_r;
            eq_A_out(idx)     = eq_iss     ./ k_cur;
            eq_dum_out(idx)   = neg;
            income_A_out(idx) = pi_all     ./ k_cur;
            firm_id_out(idx)  = ((sim-1)*N_firm+1 : sim*N_firm)';
            year_id_out(idx)  = t - (T_tot - T_keep);
            obs_ptr           = obs_ptr + N_firm;
        end

        u_z     = rand(N_firm, 1);
        iz_next = zeros(N_firm, 1);
        for iz_c = 1:Nz
            m = (iz_all == iz_c);
            if ~any(m), continue; end
            iz_next(m) = sum(bsxfun(@ge, u_z(m), Pi_cum(iz_c,:)), 2) + 1;
        end
        iz_next = min(max(iz_next, 1), Nz);

        ik_all = ikp_all;
        ip_all = ipp_all;
        iz_all = iz_next;
    end
end

n = obs_ptr;
panel.I_A      = I_A_out(1:n);
panel.EBITDA_A = EBITDA_A_out(1:n);
panel.Q        = Q_out(1:n);
panel.debt_A   = debt_A_out(1:n);
panel.eq_A     = eq_A_out(1:n);
panel.eq_dummy = eq_dum_out(1:n);
panel.income_A = income_A_out(1:n);
panel.firm_id  = firm_id_out(1:n);
panel.year_id  = year_id_out(1:n);

fprintf(['  done: %d obs | mean(I/A)=%.4f EBITDA/A=%.4f ' ...
         'debt/A=%.4f freq(eq)=%.4f Q=%.3f\n'], ...
        n, mean(panel.I_A), mean(panel.EBITDA_A), mean(panel.debt_A), ...
        mean(double(panel.eq_dummy)), mean(panel.Q));
end
