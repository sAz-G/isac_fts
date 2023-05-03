%------------------------------------------------------------------------
% FUNCTION NAME: avg_data_rate
% AUTHOR: Sharif Azem
%         Markus Krantzik
%
% DESCRIPTION: provides the average data rate for a given trajectory
%
% INPUTS:
% S      - the trajectory points
% s_c    - the position of the communication user
% params - simulation parameters and constants 
%
% OUTPUTS:
%       rt - average data rate of the given trajectory and position of the
%       communication user
%
% USAGE: grad_rate_dim = rate_grad(S,s_c,params,dim,N)
%-----------------------------------------------------------------------
function rt = avg_data_rate(S,s_c, params,N)

P              = params.sim.P;
sigma_0        = params.sim.sigma_0;
alpha_0        = params.sim.alpha_0;
B              = params.sim.B;

rel_dist       = user_quad_distance(S, s_c, params.sim.H);

snr_const = (P * alpha_0)./(sigma_0.^2);
rt = (B./N).*sum(log2(1 + snr_const./rel_dist.^2));
end

