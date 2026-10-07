function m_sim = compute_moments(panel)
% COMPUTE_MOMENTS  Simulated moments for HW2005 - revised version 2
%
%  Revision note for moments 9-10 (Serial corr): a direct AR(1) regression
%  is used instead of a differenced IV estimator. Diagnostics showed clear
%  positive autocorrelation in EBITDA_A (firm3: 0.22 -> 0.28 -> 0.29 ->
%  0.31), while the differenced IV estimator returned a negative value,
%  indicating severe small-sample bias of differenced IV on a short panel
%  with T=9. The fixed-effects demeaned direct OLS AR(1) is used instead:
%     y_dm_it = rho * y_dm_i,t-1 + u_it
%  Moment 10 (SD of shock) is estimated directly as sigma_u = std(u_it)
%  from the residuals.

m_sim = zeros(10, 1);

I_A      = panel.I_A;
EBITDA_A = panel.EBITDA_A;
Q        = panel.Q;
debt_A   = panel.debt_A;
eq_A     = panel.eq_A;
eq_dummy = double(panel.eq_dummy);
firm_id  = panel.firm_id;
year_id  = panel.year_id;

valid = isfinite(debt_A) & isfinite(Q) & isfinite(I_A) & isfinite(EBITDA_A);
if mean(~valid) > 0.01
    warning('compute_moments: %.1f%% of observations contain NaN/Inf', mean(~valid)*100);
end

%% Moments 1-6
m_sim(1) = mean(I_A(valid));
m_sim(2) = var(I_A(valid), 1);
m_sim(3) = mean(EBITDA_A(valid));
m_sim(4) = mean(debt_A(valid));
m_sim(5) = mean(eq_A(valid));
m_sim(6) = mean(eq_dummy(valid));

%% Moments 7-8: fixed-effects demeaned OLS
firms = unique(firm_id);
Nf    = length(firms);

I_dm    = zeros(size(I_A));
Q_dm    = zeros(size(Q));
debt_dm = zeros(size(debt_A));

for fi = 1:Nf
    mask = (firm_id == firms(fi)) & valid;
    if sum(mask) < 2, continue; end
    I_dm(mask)    = I_A(mask)    - mean(I_A(mask));
    Q_dm(mask)    = Q(mask)      - mean(Q(mask));
    debt_dm(mask) = debt_A(mask) - mean(debt_A(mask));
end

v2 = valid;
dQ = Q_dm(v2)' * Q_dm(v2);
if dQ > 1e-12
    m_sim(7) = (Q_dm(v2)' * I_dm(v2))    / dQ;
    m_sim(8) = (Q_dm(v2)' * debt_dm(v2)) / dQ;
else
    m_sim(7) = 0; m_sim(8) = 0;
end

%% Moments 9-10: fixed-effects demeaned direct OLS AR(1)
% y_dm_it = rho * y_dm_i,t-1 + u_it
% Using pairs of consecutive periods
income_dm = zeros(size(EBITDA_A));
for fi = 1:Nf
    mask = (firm_id == firms(fi)) & valid;
    if sum(mask) < 2, continue; end
    income_dm(mask) = EBITDA_A(mask) - mean(EBITDA_A(mask));
end

y_curr_vec = [];
y_lag_vec  = [];

for fi = 1:Nf
    mask_fi = (firm_id == firms(fi)) & valid;
    if sum(mask_fi) < 2, continue; end
    y_fi = income_dm(mask_fi);
    t_fi = year_id(mask_fi);
    [t_s, sidx] = sort(t_fi);
    y_s = y_fi(sidx);
    for tt = 2:length(t_s)
        if t_s(tt) - t_s(tt-1) == 1   % consecutive periods
            y_curr_vec(end+1,1) = y_s(tt);     %#ok<AGROW>
            y_lag_vec(end+1,1)  = y_s(tt-1);   %#ok<AGROW>
        end
    end
end

if length(y_lag_vec) > 10
    denom = y_lag_vec' * y_lag_vec;
    if denom > 1e-12
        rho_hat = (y_lag_vec' * y_curr_vec) / denom;
    else
        rho_hat = 0;
    end
    m_sim(9) = rho_hat;
    resid     = y_curr_vec - rho_hat * y_lag_vec;
    m_sim(10) = std(resid);   % direct estimate of sigma_u (not differenced, no sqrt(2) correction)
else
    warning('compute_moments: insufficient AR(1) pairs');
    m_sim(9) = 0; m_sim(10) = 0;
end

fprintf(['  moments: I/A=%.4f Var=%.5f EBITDA=%.4f debt=%.4f ' ...
         'eq=%.4f freq=%.4f IQ=%.4f DQ=%.4f AR1=%.4f SD=%.4f\n'], ...
        m_sim(1),m_sim(2),m_sim(3),m_sim(4),m_sim(5), ...
        m_sim(6),m_sim(7),m_sim(8),m_sim(9),m_sim(10));
end
