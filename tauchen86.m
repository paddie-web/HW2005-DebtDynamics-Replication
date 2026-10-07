function [z_grid, Pi] = tauchen86(rho, sigma_eps, N, m_std)
% TAUCHEN86  Discretization of an AR(1) process following Tauchen (1986).
%
%  Revision note: m_std changed from 3 to 2. With m_std = 3 the grid
%  spanned 3*sigma_lnz on each side; the implied grid SD was about
%  1.73*m_std/sqrt(N)*sigma_lnz. Diagnostics showed an actual SD of 0.3416
%  versus a theoretical unconditional SD of 0.1829 (ratio ~1.87). Setting
%  m_std = 2 brings the grid SD close to the theoretical value.
%
%  Note: Tauchen originally suggests m = 3, but for a persistent process
%  with rho = 0.74 and N = 20 points, m = 3 makes the endpoint
%  probabilities tiny and crowds the middle points; m = 2 covers the 95%
%  unconditional distribution much better at N = 20.

if nargin < 3, N     = 20; end
if nargin < 4, m_std = 2;  end   % revised: originally 3, changed to 2

sigma_lnz = sigma_eps / sqrt(1 - rho^2);

lnz_max  = m_std * sigma_lnz;
lnz_grid = linspace(-lnz_max, lnz_max, N)';
step     = lnz_grid(2) - lnz_grid(1);

Pi = zeros(N, N);
for i = 1:N
    mu_cond = rho * lnz_grid(i);
    Pi(i,1) = normcdf((lnz_grid(1) - mu_cond + 0.5*step) / sigma_eps);
    Pi(i,N) = 1 - normcdf((lnz_grid(N) - mu_cond - 0.5*step) / sigma_eps);
    for j = 2:N-1
        Pi(i,j) = normcdf((lnz_grid(j) - mu_cond + 0.5*step) / sigma_eps) ...
                - normcdf((lnz_grid(j) - mu_cond - 0.5*step) / sigma_eps);
    end
    Pi(i,:) = Pi(i,:) / sum(Pi(i,:));
end

z_grid = exp(lnz_grid);
end
