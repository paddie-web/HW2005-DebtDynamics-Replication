function grids = build_state_space(params, fixed)
% BUILD_STATE_SPACE  State space construction for HW2005 - revised version 4
%
%  Design principles for the p_grid:
%    The positive-step size must be < b*(k_mid,z_mid), so that the
%    reasonable debt region contains grid points.
%    Positive range = p_upper_max (collateral constraint upper bound)
%    Negative range = -p_upper_max (symmetric, savings side)
%    Point allocation: positive side first satisfies the step-size
%    constraint, the remainder is given to the negative side.

alpha=params(1); s_val=params(2); delta=params(3);
rho=params(4); sigma_eps=params(5);
mu_tau=params(6); sigma_tau=params(7);
r=fixed.r; Nk=fixed.Nk; Np=fixed.Np; Nz=fixed.Nz;

%% STEP 1: z grid
[z_grid, Pi] = tauchen86(rho, sigma_eps, Nz, 3);
z_max=z_grid(end); z_min=z_grid(1); z_mid=z_grid(ceil(Nz/2));

%% STEP 2: k grid
k_bar = (z_max/delta)^(1/(1-alpha));
k_ss_z1 = (alpha/delta)^(1/(1-alpha));
k_min = max(k_ss_z1/10, 1.0);
k_grid = exp(linspace(log(k_min), log(k_bar), Nk))';
k_ss_lo = (alpha*z_min/delta)^(1/(1-alpha));
k_ss_hi = (alpha*z_max/delta)^(1/(1-alpha));
fprintf('  k_grid: [%.2f,%.2f] | k_ss no-tax:[%.2f,%.2f]\n',...
        k_grid(1),k_grid(end),k_ss_lo,k_ss_hi);

%% STEP 3: collateral upper bound (face value)
p_upper = zeros(Nk,1);
for ik=1:Nk
    kv=k_grid(ik); piv=z_min*kv^alpha; yv=piv-delta*kv;
    [~,gv]=corporate_tax(yv,mu_tau,sigma_tau);
    p_upper(ik) = (s_val*kv*(1-delta)+(piv-gv))*(1+r);
end
p_upper = max(p_upper,0);
p_max_debt = max(p_upper);

%% STEP 4: compute b* (debt face value at the CF=0 boundary, approximated
%  with k_mid, z_mid, b=0)
k_mid_val = k_grid(ceil(Nk/2));
pi_mid    = z_mid * k_mid_val^alpha;
y_mid     = pi_mid - delta*k_mid_val;
[~,g_mid] = corporate_tax(y_mid, mu_tau, sigma_tau);
b_star    = (pi_mid - g_mid) * (1+r);
fprintf('  p_upper_max=%.4f | b*(k_mid,z_mid)=%.4f\n', p_max_debt, b_star);

%% STEP 5: determine positive-side step size and number of points
%  Target step = b*/3 (ensure at least 3 points in the debt region)
target_step = b_star / 3;
%  Points needed on the positive side (including the 0 endpoint)
Np_pos = ceil(p_max_debt / target_step) + 1;
%  Points on the negative side (range [-p_max_debt, 0], same step)
Np_neg = ceil(p_max_debt / target_step) + 1;

%  If the total exceeds Np, shrink proportionally (preserving step ratio)
total_needed = Np_pos + Np_neg - 1;   % -1 because 0 is shared
if total_needed > Np
    scale = (Np - 1) / (total_needed - 1);
    Np_pos = max(3, round(Np_pos * scale));
    Np_neg = Np - Np_pos + 1;
end

%  Build the piecewise grid
p_pos = linspace(0, p_max_debt, Np_pos)';
p_neg = linspace(-p_max_debt, 0, Np_neg)';
p_grid = unique([p_neg; p_pos]);   % merge and dedupe; 0 appears once

%  Trim/adjust to exactly Np points (unique may produce floating-point duplicates)
while length(p_grid) > Np
    %  Drop the negative point closest to 0 (savings side needs less resolution)
    neg_interior = find(p_grid < 0, 1, 'last');
    if ~isempty(neg_interior) && neg_interior > 1
        p_grid(neg_interior) = [];
    else
        p_grid(end-1) = [];
    end
end
while length(p_grid) < Np
    %  Insert a point at the largest positive gap
    diffs = diff(p_grid);
    [~,idx] = max(diffs);
    p_grid = sort([p_grid; mean(p_grid(idx:idx+1))]);
end

%  Ensure 0 is exactly on the grid
[~,ip_zero] = min(abs(p_grid));
p_grid(ip_zero) = 0.0;

%  Report the actual step size
pos_idx = find(p_grid > 0, 1, 'first');
actual_step_pos = p_grid(pos_idx);   % first positive step = first positive point
fprintf('  first positive step=%.4f (target < b*/3=%.4f)\n', actual_step_pos, target_step);

%% Validation
assert(all(k_grid>0),'k_grid contains non-positive values');
assert(all(diff(k_grid)>0),'k_grid is not monotonic');
assert(all(z_grid>0),'z_grid contains non-positive values');
assert(max(abs(sum(Pi,2)-1))<1e-10,'Pi rows do not sum to 1');
assert(length(p_grid)==Np, sprintf('p_grid length=%d ~= Np=%d',length(p_grid),Np));
assert(p_grid(ip_zero)==0,'p_grid(ip_zero) ~= 0');
if actual_step_pos >= b_star
    warning('first positive step %.4f >= b* %.4f; the debt region may lack grid points, consider increasing Np',...
            actual_step_pos, b_star);
end

fprintf('  p_grid:[%.2f,%.2f](%d) ip_zero=%d | p_upper:[%.2f,%.2f]\n',...
        p_grid(1),p_grid(end),Np,ip_zero,min(p_upper),max(p_upper));
fprintf('  k:[%.2f,%.2f](%d)|p:[%.2f,%.2f](%d)|z:[%.4f,%.4f](%d)\n',...
        k_grid(1),k_grid(end),Nk,p_grid(1),p_grid(end),Np,...
        z_grid(1),z_grid(end),Nz);

%% Package
grids.k_grid=k_grid; grids.p_grid=p_grid; grids.z_grid=z_grid;
grids.Pi=Pi; grids.p_upper=p_upper;
grids.Nk=Nk; grids.Np=Np; grids.Nz=Nz;
grids.k_bar=k_bar; grids.k_ss_lo=k_ss_lo; grids.k_ss_hi=k_ss_hi;
grids.p_min=p_grid(1); grids.p_max=p_grid(end); grids.ip_zero=ip_zero;
grids.alpha=alpha; grids.s_val=s_val; grids.delta=delta; grids.r=r;
end
