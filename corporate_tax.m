function [tau_c, g_y] = corporate_tax(y, mu_tau, sigma_tau)
% CORPORATE_TAX  Convex marginal tax rate and tax bill.
%
%   Paper, Assumption 6 and Section V:
%     tau_c(y) = 0.35 * Phi((y - mu_tau) / sigma_tau)
%     g(y)     = integral_0^y tau_c(x) dx  = 0.35 * [H(y) - H(0)]
%     H(x)     = (x-mu)*Phi((x-mu)/sigma) + sigma*phi((x-mu)/sigma)
%
%   Properties (satisfying Assumption 6):
%     g(0) = 0
%     lim_{y->+inf} tau_c(y) = 0.35 = tau_c_bar < 1
%     lim_{y->-inf} tau_c(y) = 0
%     tau_c_bar > tau_i  (0.35 > 0.25, bounded-savings condition holds)
%
%   For y < 0: g(y) < 0, capturing the effect of the loss carry-forward
%   provision. Supports array inputs of any shape.
%
%   INPUTS
%     y         - taxable income (scalar or array)
%     mu_tau    - location parameter of the tax-rate function (Table III: -12.267)
%     sigma_tau - scale parameter of the tax-rate function (Table III:   9.246)
%
%   OUTPUTS
%     tau_c  - marginal tax rate (same shape as y)
%     g_y    - tax bill (same shape as y)

if sigma_tau <= 0
    error('corporate_tax: sigma_tau must be positive, got %.6f', sigma_tau);
end

% Standardized variable (avoid repeated division)
u = (y - mu_tau) ./ sigma_tau;

% Marginal tax rate
tau_c = 0.35 * normcdf(u);

% Tax bill (computed only when requested)
if nargout > 1
    % H(x) = (x-mu)*Phi((x-mu)/sigma) + sigma*phi((x-mu)/sigma)
    H = @(x) (x - mu_tau) .* normcdf((x - mu_tau)./sigma_tau) ...
            + sigma_tau    .* normpdf((x - mu_tau)./sigma_tau);
    g_y = 0.35 * (H(y) - H(0));
end
end
