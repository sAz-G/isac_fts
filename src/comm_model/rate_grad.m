%------------------------------------------------------------------------
% FUNCTION NAME: rate_grad
% AUTHORS: Sharif Azem (sAz-G on GitHub) and Markus Krantzik (mardank on GitHub)
%
% DESCRIPTION: Calculates the gradient of the rate of function along a given dimension.
%
% INPUTS:
%   S         - Trajectory points at which the gradient is calculated.
%   s_c       - Position of the communication user.
%   params    - Simulation parameters and constants.
%   dim       - Dimension along which the gradient is calculated ('x' or 'y').
%   N         - Amount of trajectory points (can be different from the amount of
%               points in S, because of the multi-stage approach).
%
% OUTPUTS:
%   grad_rate_dim - Gradient in the given dimension.
%
% USAGE: grad_rate_dim = rate_grad(S, s_c, params, dim, N)
%-----------------------------------------------------------------------

function grad_rate_dim = rate_grad(S, s_c, params, dim, N)

P         = params.sim.P;
sigma_0   = params.sim.sigma_0;
alpha_0   = params.sim.alpha_0;
B         = params.sim.B;

rel_dist = user_quad_distance(S, s_c, params.sim.H);
snr_const = (P * alpha_0) / (sigma_0^2);

if strcmp(dim, 'x')
    rel_pos_x = (S(1,:) - s_c(1));
    grad_rate_dim = -2 * snr_const .* rel_pos_x ./ (rel_dist.^2 .* (rel_dist.^2 + snr_const));
    grad_rate_dim = grad_rate_dim .* B ./ (N * log(2));
elseif strcmp(dim, 'y')
    rel_pos_y = (S(2,:) - s_c(2));
    grad_rate_dim = -2 * snr_const .* rel_pos_y ./ (rel_dist.^2 .* (rel_dist.^2 + snr_const));
    grad_rate_dim = grad_rate_dim .* B ./ (N * log(2));
end
end
