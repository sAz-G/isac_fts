%------------------------------------------------------------------------
% FUNCTION NAME: rate_grad
% AUTHOR: Sharif Azem
%         Markus Krantzik
%
% DESCRIPTION:  calculates the gradient of the rate of function along a given dimension 
%
% INPUTS:
% S         - the trajectory points at which the gradient is going to be calculated
% s_c       - the position of the communication user
% params    - simulation paramseters and constants
% dim       - x or y dimension.
% N         -  amount of trajectory points (can be different from the amount of
%               points in S, because of the multi stage approach).
%
% OUTPUTS:
%       grad_rate_dim - gradient in the given dimension
%
% USAGE: grad_rate_dim = rate_grad(S,s_c,params,dim,N)
%-----------------------------------------------------------------------

function grad_rate_dim = rate_grad(S,s_c,params,dim,N)

P              = params.sim.P;
sigma_0        = params.sim.sigma_0;
alpha_0        = params.sim.alpha_0;
B              = params.sim.B;

rel_dist       = user_quad_distance(S, s_c, params.sim.H);
snr_const      = (P.*alpha_0)./(sigma_0.^2);
if strcmp(dim,'x')
    rel_pos_x     = (S(1,:) - s_c(1));
    grad_rate_dim = -2*snr_const.*rel_pos_x./( rel_dist.^2.*( rel_dist.^2 + snr_const));
    grad_rate_dim = grad_rate_dim.*B./(N * log(2));
elseif strcmp(dim,'y')
    rel_pos_y     = (S(2,:) - s_c(2));
    grad_rate_dim = -2*snr_const.*rel_pos_y./( rel_dist.^2.*( rel_dist.^2 + snr_const));
    grad_rate_dim = grad_rate_dim.*B./(N * log(2));
end
end