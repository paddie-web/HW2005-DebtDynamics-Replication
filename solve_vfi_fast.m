function [V, k_idx_pol, p_idx_pol, vfi_history] = solve_vfi_fast(grids, params, fixed, vfi_opts)
% SOLVE_VFI_FAST  Value function iteration for HW2005.
%  Returns an additional vfi_history struct for the reporting script:
%    vfi_history.diff_vec  : sup-norm convergence value in each iteration
%    vfi_history.iter_final: final number of iterations
%    vfi_history.converged : whether VFI converged
%    vfi_history.tol       : convergence threshold
%    vfi_history.elapsed_s : total elapsed time (seconds)

alpha=params(1); s_val=params(2); delta=params(3);
mu_tau=params(6); sigma_tau=params(7); lambda=params(8);
r=fixed.r; tau_d=fixed.tau_d; tau_i=fixed.tau_i;
beta = 1/(1+r*(1-tau_i));

k_grid=grids.k_grid; p_grid=grids.p_grid; z_grid=grids.z_grid;
Pi=grids.Pi; p_upper=grids.p_upper;
Nk=grids.Nk; Np=grids.Np; Nz=grids.Nz;

if nargin<4||isempty(vfi_opts)
    vfi_opts.tol=1e-6; vfi_opts.max_iter=2000;
end
tol=vfi_opts.tol; max_iter=vfi_opts.max_iter;
inv1r = 1/(1+r);

k3=reshape(k_grid,[Nk,1,1]);
p3=reshape(p_grid,[1,Np,1]);
z3=reshape(z_grid,[1,1,Nz]);

pi_kz   = bsxfun(@times, z3, bsxfun(@power,k3,alpha));
pi_full = repmat(pi_kz,[1,Np,1]);
y_kpz   = pi_full - delta*k3 - r*p3*inv1r;
[~,g_kpz]   = corporate_tax(y_kpz, mu_tau, sigma_tau);
pi_net_kpz  = pi_full - g_kpz;
bv3         = p3 * inv1r;
CF_kpz      = pi_net_kpz - bv3;
Phi_s       = double(pi_net_kpz < bv3);
phi_all     = (1-Phi_s) + s_val*Phi_s;
n_kpz       = Phi_s .* max(0, (bv3-pi_net_kpz)/s_val);
k_eff_kpz   = k3*(1-delta) - n_kpz;
term1_kpz   = phi_all .* CF_kpz;
pp_3d       = reshape(p_grid*inv1r, [1,1,Np]);
feasible    = bsxfun(@le, reshape(p_grid,[1,Np]), reshape(p_upper,[Nk,1])+1e-10);

fprintf('VFI: %dx%dx%d  beta=%.6f  tol=%.0e  max_iter=%d\n',...
        Nk,Np,Nz,beta,tol,max_iter);

V_init = pi_kz * (1-tau_d) / (1-beta);
V = repmat(V_init,[1,Np,1]);
k_idx_pol = ones(Nk,Np,Nz,'int32');
p_idx_pol = ones(Nk,Np,Nz,'int32');
diff_val  = Inf;

% Record the convergence history
diff_vec = zeros(max_iter,1);

t_start = tic;
for iter = 1:max_iter
    EV = reshape(reshape(V,Nk*Np,Nz)*Pi', Nk,Np,Nz);
    V_new     = -Inf(Nk,Np,Nz);
    k_pol_new = ones(Nk,Np,Nz,'int32');
    p_pol_new = ones(Nk,Np,Nz,'int32');

    for iz = 1:Nz
        t1_iz   = term1_kpz(:,:,iz);
        keff_iz = k_eff_kpz(:,:,iz);
        EV_iz   = EV(:,:,iz);
        for ikp = 1:Nk
            kp_val  = k_grid(ikp);
            inv_mat = kp_val - keff_iz;
            t1_3d   = reshape(t1_iz - inv_mat,[Nk,Np,1]);
            bracket_3d = bsxfun(@plus, t1_3d, pp_3d);
            e_3d    = zeros(Nk,Np,Np);
            pos3    = bracket_3d >  0;
            neg3    = bracket_3d <  0;
            e_3d(pos3) = (1-tau_d)  .* bracket_3d(pos3);
            e_3d(neg3) = (1+lambda) .* bracket_3d(neg3);
            EV_row  = reshape(EV_iz(ikp,:),[1,1,Np]);
            rhs_3d  = e_3d + beta*EV_row;
            feas_row = reshape(feasible(ikp,:),[1,1,Np]);
            rhs_3d(~repmat(feas_row,[Nk,Np,1])) = -Inf;
            [best_val, best_ipp] = max(rhs_3d,[],3);
            improve = best_val > V_new(:,:,iz);
            if any(improve(:))
                tmp_V=V_new(:,:,iz); tmp_k=k_pol_new(:,:,iz); tmp_p=p_pol_new(:,:,iz);
                tmp_V(improve)=best_val(improve);
                tmp_k(improve)=int32(ikp);
                tmp_p(improve)=int32(best_ipp(improve));
                V_new(:,:,iz)=tmp_V; k_pol_new(:,:,iz)=tmp_k; p_pol_new(:,:,iz)=tmp_p;
            end
        end
    end

    diff_val = max(abs(V_new(:)-V(:)));
    diff_vec(iter) = diff_val;   % store the diff of this iteration
    V=V_new; k_idx_pol=k_pol_new; p_idx_pol=p_pol_new;

    if mod(iter,50)==0||diff_val<tol*10
        fprintf('  iter=%4d | diff=%.3e | %.1fs\n',iter,diff_val,toc(t_start));
    end
    if diff_val<tol
        fprintf('  VFI converged: iter=%d diff=%.3e %.1fs\n',iter,diff_val,toc(t_start));
        diff_vec = diff_vec(1:iter);
        break;
    end
end
if diff_val>=tol
    warning('VFI did not converge: diff=%.3e',diff_val);
    diff_vec = diff_vec(1:max_iter);
end

elapsed = toc(t_start);
fprintf('  V range: [%.3f,%.3f] | k_pol:[%d,%d] | p_pol:[%d,%d]\n',...
        min(V(:)),max(V(:)),min(k_idx_pol(:)),max(k_idx_pol(:)),...
        min(p_idx_pol(:)),max(p_idx_pol(:)));
fprintf('  fraction p_pol>0 (debt): %.1f%%\n', mean(p_idx_pol(:)>grids.ip_zero)*100);

% Package the history
vfi_history.diff_vec   = diff_vec;
vfi_history.iter_final = length(diff_vec);
vfi_history.converged  = (diff_val < tol);
vfi_history.tol        = tol;
vfi_history.elapsed_s  = elapsed;
vfi_history.beta       = beta;
vfi_history.Nk = Nk; vfi_history.Np = Np; vfi_history.Nz = Nz;
end
