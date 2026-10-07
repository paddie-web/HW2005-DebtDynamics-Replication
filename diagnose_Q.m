% DIAGNOSE_Q.M - Diagnostic script for the Investment-Q sensitivity sign
% reversal. Requires workspace variables: panel, grids, V, k_idx_pol,
% p_idx_pol, b_paper, fixed.

fprintf('\n===== A: Basic Q statistics =====\n');
fprintf('mean(Q)   = %.4f  (paper approx 1-2)\n', mean(panel.Q));
fprintf('median(Q) = %.4f\n', median(panel.Q));
fprintf('std(Q)    = %.4f\n', std(panel.Q));
fprintf('min(Q)    = %.4f\n', min(panel.Q));
fprintf('max(Q)    = %.4f\n', max(panel.Q));
fprintf('fraction Q<0  = %.1f%%\n', mean(panel.Q<0)*100);
fprintf('fraction Q>10 = %.1f%%\n', mean(panel.Q>10)*100);

fprintf('\n===== B: Raw correlations of Q with investment =====\n');
valid = isfinite(panel.Q) & isfinite(panel.I_A);
fprintf('corr(Q, I/A)        = %.4f\n', corr(panel.Q(valid), panel.I_A(valid)));
fprintf('corr(Q, EBITDA/A)   = %.4f\n', corr(panel.Q(valid), panel.EBITDA_A(valid)));
fprintf('corr(Q, debt/A)     = %.4f\n', corr(panel.Q(valid), panel.debt_A(valid)));

fprintf('\n===== C: Correlations after demeaning by firm =====\n');
firms = unique(panel.firm_id);
Nf = length(firms);
I_dm = zeros(size(panel.I_A));
Q_dm = zeros(size(panel.Q));
for fi = 1:Nf
    mask = panel.firm_id==firms(fi) & valid;
    if sum(mask)<2, continue; end
    I_dm(mask) = panel.I_A(mask) - mean(panel.I_A(mask));
    Q_dm(mask) = panel.Q(mask)   - mean(panel.Q(mask));
end
fprintf('corr(Q_dm, I_dm) = %.4f  (Investment-Q sensitivity should be positive)\n',...
    corr(Q_dm(valid), I_dm(valid)));

fprintf('\n===== D: Decomposition of Q (a few sample states) =====\n');
fprintf('Q = (EV(k'',b'',z) + b''/(1+r)) / k''\n');
fprintf('EV = E_z[V(k'',b'',z'')|z]\n\n');

r=fixed.r; alpha=b_paper(1);
iz_mid=ceil(grids.Nz/2);
ik_mid=ceil(grids.Nk/2);
ip0=grids.ip_zero;

V_mat = reshape(V, grids.Nk*grids.Np, grids.Nz);
EV_mat = V_mat * grids.Pi';

fprintf('fix (k_mid,b=0), scan z:\n');
fprintf('  iz | z_val | k_next | b_next | EV     | bv_next | Q\n');
for iz=1:grids.Nz
    lin=sub2ind([grids.Nk,grids.Np,grids.Nz],ik_mid,ip0,iz);
    ikp=double(k_idx_pol(lin));
    ipp=double(p_idx_pol(lin));
    k_next=grids.k_grid(ikp);
    b_next=grids.p_grid(ipp);
    bv_next=b_next/(1+r);
    lin_kp=(ikp-1)*grids.Np+ipp;
    ev=EV_mat(lin_kp, iz);
    Q_val=(ev+bv_next)/k_next;
    fprintf('  %2d | %.4f | %.4f  | %.4f  | %.4f | %.4f   | %.4f\n',...
        iz,grids.z_grid(iz),k_next,b_next,ev,bv_next,Q_val);
end

fprintf('\n===== E: What Q should look like theoretically =====\n');
fprintf('Paper definition: Q_t = (market value of equity + market value of debt) / capital\n');
fprintf('V(k,b,z) is equity value (discounted stream of future dividends)\n');
fprintf('Q = (V + b/(1+r)) / k\n');
fprintf('When z is high, expected future profitability is high -> V high -> Q high -> more investment\n');
fprintf('Therefore corr(Q_dm, I_dm) should be positive\n');
fprintf('\nCurrent pattern of Q vs z (from the table above):\n');
fprintf('If Q is monotonically increasing in z, the Q calculation is correct and the problem lies elsewhere\n');
fprintf('If Q is non-monotonic or decreasing in z, the Q calculation is wrong\n');

fprintf('\nDiagnostics complete.\n');
