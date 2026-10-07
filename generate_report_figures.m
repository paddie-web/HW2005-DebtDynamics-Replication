%% GENERATE_REPORT_FIGURES.M
%  Run after run_hw2005_with_paper_params.m completes.
%  Requires workspace: V, k_idx_pol, p_idx_pol, vfi_history,
%                      grids, panel, m_sim, M_hat, b_paper, fixed

close all;

c1=[0.122 0.467 0.706]; c2=[0.839 0.153 0.157];
c3=[0.173 0.627 0.173]; c4=[0.580 0.404 0.741];
iz_lo=3; iz_mid=ceil(grids.Nz/2); iz_hi=grids.Nz-2;
out_prefix='report';

fprintf('============================================================\n');
fprintf('  Generating report figures (HW2005 replication)\n');
fprintf('============================================================\n\n');

%% ================================================================
%  Fig1: VFI convergence curve
%  Revision: ylabel uses plain text instead of LaTeX to avoid
%  interpreter errors
%% ================================================================
fig1=figure('Name','Fig1_VFI_Convergence','Position',[50,550,700,380]);
semilogy(1:vfi_history.iter_final, vfi_history.diff_vec,...
         'Color',c1,'LineWidth',1.6);
hold on;
yline(vfi_history.tol,'--r','LineWidth',1.2);
text(vfi_history.iter_final*0.6, vfi_history.tol*3,...
     sprintf('tol = %.0e',vfi_history.tol),'Color',c2,'FontSize',9);

xlabel('Iteration','FontSize',11);
ylabel('Sup-norm error ||V_{n+1} - V_n||','FontSize',11);   % plain text, no LaTeX
title(sprintf('VFI convergence  (%d iterations,  %.0f s)',...
      vfi_history.iter_final, vfi_history.elapsed_s),'FontSize',12);

conv_str = sprintf('final error: %.2e\nbeta = %.4f\nstate space: %dx%dx%d',...
    vfi_history.diff_vec(end), vfi_history.beta,...
    vfi_history.Nk, vfi_history.Np, vfi_history.Nz);
text(vfi_history.iter_final*0.52, vfi_history.diff_vec(1)*0.25, conv_str,...
     'FontSize',9,'BackgroundColor','w','EdgeColor',[0.7 0.7 0.7]);
grid on; box on;
saveas(fig1,[out_prefix '_fig1_convergence.png']);
fprintf('[Fig1] VFI convergence curve saved\n');

%% ================================================================
%  Fig2: Stationary distribution
%% ================================================================
fig2=figure('Name','Fig2_Distribution','Position',[50,120,900,380]);

k_pol_dist=double(k_idx_pol(:));
p_pol_dist=double(p_idx_pol(:));

subplot(1,2,1);
k_vals=grids.k_grid(k_pol_dist);
histogram(k_vals,30,'Normalization','probability',...
          'FaceColor',c1,'EdgeColor','none','FaceAlpha',0.75);
xlabel('Next-period capital k''','FontSize',10);
ylabel('Probability','FontSize',10);
title('Capital policy distribution k''(k,p,z)','FontSize',11);
hold on;
xline(grids.k_ss_lo,'--','Color',c2,'LineWidth',1.2,'Label','k_{ss}(z_{min})');
xline(grids.k_ss_hi,'--','Color',c3,'LineWidth',1.2,'Label','k_{ss}(z_{max})');
grid on; box on;

subplot(1,2,2);
p_vals=grids.p_grid(p_pol_dist);
histogram(p_vals,30,'Normalization','probability',...
          'FaceColor',c2,'EdgeColor','none','FaceAlpha',0.75);
xlabel('Next-period debt face value b''','FontSize',10);
ylabel('Probability','FontSize',10);
title('Debt policy distribution b''(k,p,z)','FontSize',11);
hold on;
xline(0,'k-','LineWidth',1.5,'Label','b=0');
grid on; box on;

sgtitle('Stationary distribution of policy functions (across all states)','FontSize',12,'FontWeight','bold');
saveas(fig2,[out_prefix '_fig2_distribution.png']);
fprintf('[Fig2] stationary distribution figure saved\n');

%% ================================================================
%  Fig3: Summary table of key steady-state variables
%  Revision: annotation textbox used instead of line + Units=normalized
%% ================================================================
fig3=figure('Name','Fig3_Summary_Table','Position',[760,550,640,360]);
axes('Position',[0 0 1 1],'Visible','off');
hold on;

mean_Q      = mean(panel.Q);
mean_IA     = mean(panel.I_A);
mean_EBITDA = mean(panel.EBITDA_A);
mean_debt   = mean(panel.debt_A);
freq_eq     = mean(double(panel.eq_dummy));
pct_save    = mean(panel.debt_A<0)*100;
beta_val    = vfi_history.beta;
Vmin=min(V(:)); Vmax=max(V(:));
k_pol_min=grids.k_grid(min(k_idx_pol(:)));
k_pol_max=grids.k_grid(max(k_idx_pol(:)));

col1={'Solution method';'Discount factor beta';'Iterations to convergence';'Final error';...
      'V range';'k policy range';...
      'Mean I/A';'Mean EBITDA/A';'Mean net debt/A';...
      'Equity issuance frequency';'Fraction of saving firms';'Mean Tobin Q'};
col2={
    'Howard VFI (no acceleration)';
    sprintf('%.6f', beta_val);
    sprintf('%d iters', vfi_history.iter_final);
    sprintf('%.2e', vfi_history.diff_vec(end));
    sprintf('[%.1f, %.1f]', Vmin, Vmax);
    sprintf('[%.1f, %.1f]', k_pol_min, k_pol_max);
    sprintf('%.4f  (paper: 0.079)', mean_IA);
    sprintf('%.4f  (paper: 0.146)', mean_EBITDA);
    sprintf('%.4f  (paper: 0.075)', mean_debt);
    sprintf('%.4f  (paper: 0.099)', freq_eq);
    sprintf('%.1f%%  (paper approx 10%%)', pct_save);
    sprintf('%.3f  (paper approx 1-2)', mean_Q);
};

nrow=length(col1);
% Title
text(0.5,0.96,'HW2005 key steady-state variables','FontSize',13,'FontWeight','bold',...
     'HorizontalAlignment','center','Units','normalized');
% Separator lines via annotation
annotation('line',[0.02 0.98],[0.91 0.91],'Color',[0.3 0.3 0.3],'LineWidth',1.2);
annotation('line',[0.02 0.98],[0.03 0.03],'Color',[0.3 0.3 0.3],'LineWidth',0.8);
annotation('line',[0.40 0.40],[0.03 0.91],'Color',[0.8 0.8 0.8],'LineWidth',0.8);

for rr=1:nrow
    yp=0.88 - (rr-1)*(0.85/nrow);
    text(0.03, yp, col1{rr},'FontSize',9,'FontWeight','bold',...
         'Units','normalized','HorizontalAlignment','left');
    text(0.42, yp, col2{rr},'FontSize',9,...
         'Units','normalized','HorizontalAlignment','left');
end
saveas(fig3,[out_prefix '_fig3_summary_table.png']);
fprintf('[Fig3] steady-state summary table saved\n');

%% ================================================================
%  Fig4: Policy function k'(k, z)
%% ================================================================
fig4=figure('Name','Fig4_Policy_Capital','Position',[760,120,700,400]);
ip_fix=grids.ip_zero;
iz_labels={'low z','mid z','high z'};
iz_vals=[iz_lo, iz_mid, iz_hi];
colors_4={c1,c4,c2};

hold on;
for jj=1:3
    iz=iz_vals(jj);
    kp_idx=double(squeeze(k_idx_pol(:,ip_fix,iz)));
    kp_val=grids.k_grid(kp_idx);
    plot(grids.k_grid, kp_val,'Color',colors_4{jj},'LineWidth',1.8,...
         'DisplayName',iz_labels{jj});
end
plot(grids.k_grid,grids.k_grid,'k--','LineWidth',0.9,'DisplayName','45-degree line');

xlabel('Current capital k','FontSize',11);
ylabel('Next-period capital k''','FontSize',11);
title('Capital investment policy function k''(k, b=0, z)','FontSize',12);
legend('Location','northwest','FontSize',10);
xlim([grids.k_grid(1), grids.k_grid(end)]);
grid on; box on;
saveas(fig4,[out_prefix '_fig4_policy_capital.png']);
fprintf('[Fig4] capital policy function figure saved\n');

%% ================================================================
%  Fig5: Policy function b'(k, z)
%% ================================================================
fig5=figure('Name','Fig5_Policy_Debt','Position',[50,120,700,400]);
hold on;
for jj=1:3
    iz=iz_vals(jj);
    pp_idx=double(squeeze(p_idx_pol(:,ip_fix,iz)));
    pp_val=grids.p_grid(pp_idx);
    plot(grids.k_grid, pp_val,'Color',colors_4{jj},'LineWidth',1.8,...
         'DisplayName',iz_labels{jj});
end
yline(0,'k--','LineWidth',1.2,'Label','b=0');

xlabel('Current capital k','FontSize',11);
ylabel('Next-period debt face value b''','FontSize',11);
title('Debt policy function b''(k, b=0, z)','FontSize',12);
legend('Location','best','FontSize',10);
xlim([grids.k_grid(1), grids.k_grid(end)]);
grid on; box on;
saveas(fig5,[out_prefix '_fig5_policy_debt.png']);
fprintf('[Fig5] debt policy function figure saved\n');

%% ================================================================
%  Fig6: Table II simulated-moments comparison
%  Revision: annotation used instead of line, Units=normalized removed
%% ================================================================
fig6=figure('Name','Fig6_Moments_Table','Position',[50,120,940,430]);
axes('Position',[0 0 1 1],'Visible','off');
hold on;

% Simulated moments column of Table II of the paper (values as printed)
M_paper_sim=[0.107;0.014;0.165;0.127;0.052;0.091;...
             0.036;-0.219;0.598;0.120];
moment_names={'Mean(I/A)';'Var(I/A)';'Mean(EBITDA/A)';...
              'Mean(net debt/A)';'Mean(eq issuance/A)';'Freq(eq issuance)';...
              'Investment-Q sens';'Debt-Q sens';...
              'Serial corr(income/A)';'SD(shock income/A)'};

hdrs={'Moment','Data (Table II)','Paper simulated','This replication','Diff% (vs data)'};
hx=[0.01,0.30,0.47,0.62,0.79];

text(0.5,0.97,'Table II simulated-moments comparison (HW2005)','FontSize',12,...
     'FontWeight','bold','HorizontalAlignment','center','Units','normalized');
for hh=1:length(hdrs)
    text(hx(hh),0.91,hdrs{hh},'FontSize',9,'FontWeight','bold',...
         'Units','normalized');
end
annotation('line',[0.01 0.99],[0.88 0.88],'Color',[0.3 0.3 0.3],'LineWidth',1.2);

for mm=1:10
    yp=0.85-(mm-1)*0.079;
    act=M_hat(mm); psim=M_paper_sim(mm); mysim=m_sim(mm);
    if abs(act)>1e-8
        pct=(mysim-act)/abs(act)*100;
    else
        pct=NaN;
    end
    if ~isnan(pct)&&abs(pct)>20,  fc=[0.85 0.1 0.1];
    elseif ~isnan(pct)&&abs(pct)>5, fc=[0.85 0.55 0.0];
    else,  fc=[0.1 0.5 0.1]; end

    text(hx(1),yp,moment_names{mm},'FontSize',8.5,'Units','normalized');
    text(hx(2),yp,sprintf('%.4f',act),'FontSize',8.5,'Units','normalized');
    text(hx(3),yp,sprintf('%.4f',psim),'FontSize',8.5,'Units','normalized');
    text(hx(4),yp,sprintf('%.4f',mysim),'FontSize',8.5,'Units','normalized','Color',fc);
    if isnan(pct)
        text(hx(5),yp,'N/A','FontSize',8.5,'Units','normalized','Color',fc);
    else
        text(hx(5),yp,sprintf('%.1f%%',pct),'FontSize',8.5,'Units','normalized','Color',fc);
    end
end

annotation('line',[0.01 0.99],[0.06 0.06],'Color',[0.3 0.3 0.3],'LineWidth',0.8);
obj_val=(M_hat-m_sim)'*diag(1./(M_hat.^2+1e-8))*(M_hat-m_sim);
max_pct=max(abs((m_sim-M_hat)./(abs(M_hat)+1e-8)*100));
text(0.5,0.03,sprintf('GMM objective value = %.4f   max relative deviation = %.1f%%',...
     obj_val, max_pct),'FontSize',8.5,'Units','normalized',...
     'HorizontalAlignment','center','Color',[0.3 0.3 0.3]);

saveas(fig6,[out_prefix '_fig6_moments_table.png']);
fprintf('[Fig6] moments comparison table saved\n');

%% ================================================================
fprintf('\n============================================================\n');
fprintf('  All figures generated (prefix: %s_figX_*.png)\n', out_prefix);
fprintf('  Fig1: VFI convergence curve\n');
fprintf('  Fig2: policy distributions (k and b)\n');
fprintf('  Fig3: key steady-state variables summary table\n');
fprintf('  Fig4: capital policy function k(k,z)\n');
fprintf('  Fig5: debt policy function b(k,z)\n');
fprintf('  Fig6: Table II simulated-moments comparison\n');
fprintf('============================================================\n');
